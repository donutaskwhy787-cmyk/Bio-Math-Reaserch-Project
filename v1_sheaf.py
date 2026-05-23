"""
PINGO v1 — corrected sheaf engine + fixed optimal-transport robustness test.

Fixes vs. the original Main.py prototype:
  (1) Data loading: name-based columns, real type casting, correct categorical
      encodings (IDH 'Mutant/WT', MGMT 'Methylated/Unmethylated', EGFR 'amp/normal',
      grade G2/G3/G4 -> 2/3/4, Status '1:DECEASED'->1). Missingness imputed + flagged.
  (2) A1 leakage fix: survival (Months, Status) is HELD OUT of the sheaf stalks.
      It is never used to build SRIS — only kept aside as the future endpoint.
  (3) Sheaf: 3-node modality sheaf (D=genomic, R=regulatory, C=baseline-clinical),
      identity-sheaf v1 with 1-D standardized "aggressiveness" stalks.
      SRIS(p) = ||r_DR||^2 + ||r_DC||^2 + ||r_RC||^2 = x_p^T L_F x_p.
  (4) Optimal transport: REPLACES the broken ot.solve_sample(c, c) self-coupling.
      Now a subtype->subtype transport (IDH-mut -> IDH-wt) via entropic OT (Sinkhorn),
      reporting the Wasserstein cost and per-edge residual stability under that shift.

Runs on numpy + scipy only.
"""
import numpy as np
from scipy import stats

NA = {"NA", "None", "nan", "", "NaN"}

# ---------------------------------------------------------------- load & clean
def load(path="data.txt"):
    with open(path) as f:
        header = f.readline().split()
        rows = [ln.split() for ln in f if ln.strip()]
    idx = {h: i for i, h in enumerate(header)}
    def col(name):
        i = idx[name]
        return [r[i] if i < len(r) else "NA" for r in rows]
    return col, len(rows)

def as_float(vals):
    out = np.array([np.nan if v in NA else float(v) for v in vals], float)
    return out

def as_binary(vals, one_label, zero_label):
    out = []
    for v in vals:
        if v in NA:            out.append(np.nan)
        elif v == one_label:   out.append(1.0)
        elif v == zero_label:  out.append(0.0)
        else:                  out.append(np.nan)   # unexpected token -> missing
    return np.array(out, float)

def impute_z(x, sign=+1.0):
    """median-impute, then z-score, then apply sign (so + = more aggressive)."""
    x = x.copy()
    miss = np.isnan(x)
    if miss.any():
        x[miss] = np.nanmedian(x)
    sd = x.std()
    z = (x - x.mean()) / (sd if sd > 0 else 1.0)
    return sign * z, miss

# ---------------------------------------------------------------- main
col, n = load()
print(f"Loaded {n} patients\n")

# ---- survival HELD OUT (A1): parsed but never enters the sheaf -------------
os_months = as_float(col("Months"))
deceased  = as_binary(col("Status"), "1:DECEASED", "0:LIVING")

# ---- grouping labels (for validation / OT, NOT stalk inputs) ---------------
idh_label   = col("IDH")                    # Mutant / WT
codel_label = col("IDH/codelsubtype")       # IDHmut-codel / IDHmut-non-codel / IDHwt
grade_int   = as_float([{"G2":"2","G3":"3","G4":"4"}.get(g, "NA")
                        for g in col("NeoplasmHistologicGrade")])

# ---- NODE D: genomic (sign = direction of MORE aggressive biology) ---------
idh_b   = as_binary(col("IDH"), "Mutant", "WT")          # mutant = better -> sign -
mgmt_b  = as_binary(col("MGMT"), "Methylated", "Unmethylated")  # methyl = better -> -
atrx_b  = as_binary(col("ATRXstatus"), "Mutant", "WT")   # ATRXmut ~ better -> -
mutc    = as_float(col("MutationCount"))
aneu    = as_float(col("Percentaneuploidy"))
chr710  = as_binary(col("Chr7gain/Chr10loss"), "Gainchr7&losschr10", "NocombinedCNA")  # GBM hallmark +

D_feats = [impute_z(idh_b, -1)[0], impute_z(mgmt_b, -1)[0], impute_z(atrx_b, -1)[0],
           impute_z(mutc, +1)[0],  impute_z(aneu, +1)[0],   impute_z(chr710, +1)[0]]

# ---- NODE R: regulatory / transcriptomic ----------------------------------
egfr_b  = as_binary(col("EGFR"), "amp", "normal")        # amp = aggressive +
tert_e  = as_float(col("TERTexpression(log2)"))          # high TERT = aggressive +
imm     = as_float(col("ESTIMATEimmunescore"))           # higher infiltration ~ mesenchymal +
strm    = as_float(col("ESTIMATEstromalscore"))          # +

R_feats = [impute_z(egfr_b, +1)[0], impute_z(tert_e, +1)[0],
           impute_z(imm, +1)[0],    impute_z(strm, +1)[0]]

# ---- NODE C: baseline clinical (NO survival — A1) -------------------------
age     = as_float(col("DiagnosisAge"))                  # older = worse +
kps     = as_float(col("KarnofskyPerformanceScore"))     # higher KPS = better -> -
C_feats = [impute_z(grade_int, +1)[0], impute_z(age, +1)[0], impute_z(kps, -1)[0]]

# ---- 1-D node "aggressiveness" stalks, then standardize so D,R,C compare ----
def composite(feat_list):
    s = np.mean(np.vstack(feat_list), axis=0)
    return (s - s.mean()) / (s.std() if s.std() > 0 else 1.0)

xD, xR, xC = composite(D_feats), composite(R_feats), composite(C_feats)

# ---- sheaf residuals + SRIS (identity sheaf over the D-R-C triangle) -------
r_DR = xD - xR
r_DC = xD - xC
r_RC = xR - xC
SRIS = r_DR**2 + r_DC**2 + r_RC**2
eps = 1e-9
phi_DR = r_DR**2 / (SRIS + eps)
phi_DC = r_DC**2 / (SRIS + eps)
phi_RC = r_RC**2 / (SRIS + eps)

# Equivalent Laplacian form check: L_F for a triangle (1-D stalks, identity maps)
L_F = np.array([[2.,-1.,-1.],[-1.,2.,-1.],[-1.,-1.,2.]])
X = np.vstack([xD, xR, xC]).T                      # (n,3)
SRIS_lap = np.einsum("bi,ij,bj->b", X, L_F, X) / 1.0
# triangle coboundary gives ||Bx||^2 = sum of squared edge diffs == SRIS
assert np.allclose(SRIS, np.einsum("bi,ij,bj->b", X, L_F, X)), "Laplacian/SRIS mismatch"
print("Sheaf identity check passed:  SRIS == x^T L_F x\n")

print("SRIS distribution: "
      f"median={np.median(SRIS):.3f}  mean={SRIS.mean():.3f}  max={SRIS.max():.3f}")
print(f"Mean edge fractions:  DR={phi_DR.mean():.2f}  DC={phi_DC.mean():.2f}  RC={phi_RC.mean():.2f}\n")

# ---------------------------------------------------------------- validation
def mwu(mask_a, mask_b, label):
    a, b = SRIS[mask_a], SRIS[mask_b]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    print(f"  {label:28s} median {np.median(a):.3f} vs {np.median(b):.3f}   p={p:.2e}  (n={mask_a.sum()},{mask_b.sum()})")

idh_arr = np.array(idh_label)
print("SRIS group comparisons (A1: survival not in score):")
mwu(idh_arr == "Mutant", idh_arr == "WT", "IDH mutant vs wildtype")

g = grade_int
gg = [SRIS[g == k] for k in (2, 3, 4) if (g == k).any()]
H, p = stats.kruskal(*gg)
print(f"  {'SRIS by grade G2/G3/G4':28s} medians "
      + "/".join(f"{np.median(SRIS[g==k]):.2f}" for k in (2,3,4)) + f"   KW p={p:.2e}")

cl = np.array(codel_label)
groups = ["IDHmut-codel", "IDHmut-non-codel", "IDHwt"]
gg = [SRIS[cl == name] for name in groups if (cl == name).any()]
H, p = stats.kruskal(*gg)
print(f"  {'SRIS by IDH/codel subtype':28s} "
      + "/".join(f"{name}:{np.median(SRIS[cl==name]):.2f}" for name in groups) + f"   KW p={p:.2e}\n")

# ---------------------------------------------------------------- FIXED OT
def sinkhorn(M, reg=0.05, n_iter=1000):
    """Entropic OT with uniform marginals. Returns transport plan P."""
    M = M / (M.max() + 1e-12)
    K = np.exp(-M / reg)
    a = np.full(M.shape[0], 1.0 / M.shape[0])
    b = np.full(M.shape[1], 1.0 / M.shape[1])
    u, v = np.ones_like(a), np.ones_like(b)
    for _ in range(n_iter):
        u = a / (K @ v + 1e-12)
        v = b / (K.T @ u + 1e-12)
    return u[:, None] * K * v[None, :]

# molecular feature space (D+R inputs only; no survival) for the cost matrix
Z = np.vstack(D_feats + R_feats).T                 # (n, features)
A_mask = idh_arr == "Mutant"
B_mask = idh_arr == "WT"
ZA, ZB = Z[A_mask], Z[B_mask]

# cost = squared Euclidean distance between an IDH-mut and an IDH-wt patient
M = ((ZA[:, None, :] - ZB[None, :, :]) ** 2).sum(-1)
P = sinkhorn(M)
wass = float((P * M).sum())                        # transported cost (Wasserstein-like)

# baseline: cost of a random (uniform) coupling, for reference
rand_cost = float(M.mean())
print("FIXED optimal transport  (IDH-mutant  ->  IDH-wildtype):")
print(f"  groups: {A_mask.sum()} mutant  ->  {B_mask.sum()} wildtype")
print(f"  transported (Wasserstein) cost = {wass:.3f}   vs random coupling = {rand_cost:.3f}")
print(f"  -> OT finds matches {100*(1-wass/rand_cost):.0f}% cheaper than random "
      f"(real cross-subtype structure, not the old identity self-map)\n")

# residual stability under the subtype shift: do the edge residuals line up
# between an IDH-mut patient and the IDH-wt patient OT maps it to?
def ot_stability(res):
    rA, rB = res[A_mask], res[B_mask]
    diff = np.abs(rA[:, None] - rB[None, :])
    return 1.0 - float((P * diff).sum() / (P.sum() + 1e-12)) / (np.abs(res).mean() + 1e-9)

print("OT residual stability under IDH-mut -> IDH-wt shift (1 = perfectly preserved):")
for name, res in [("r_DR", r_DR), ("r_DC", r_DC), ("r_RC", r_RC)]:
    print(f"  {name}: {ot_stability(res):+.3f}")

# ---------------------------------------------------------------- write out
import os
os.makedirs("results", exist_ok=True)
with open("results/v1_sris.csv", "w") as f:
    f.write("patient_id,SRIS,r_DR,r_DC,r_RC,phi_DR,phi_DC,phi_RC,idh,codel,grade,os_months_HELDOUT,deceased_HELDOUT\n")
    pid = col("PatientID")
    for i in range(n):
        f.write(f"{pid[i]},{SRIS[i]:.4f},{r_DR[i]:.4f},{r_DC[i]:.4f},{r_RC[i]:.4f},"
                f"{phi_DR[i]:.3f},{phi_DC[i]:.3f},{phi_RC[i]:.3f},{idh_label[i]},"
                f"{codel_label[i]},{grade_int[i]:.0f},{os_months[i]:.1f},{deceased[i]:.0f}\n")
np.savetxt("results/v1_ot_plan_mut_to_wt.txt", P, fmt="%.6e")
print("\nWrote results/v1_sris.csv and results/v1_ot_plan_mut_to_wt.txt")
