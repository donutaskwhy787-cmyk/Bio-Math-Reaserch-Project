# Proposed Additions to the BIBM Sheaf Blueprint

These three amendments keep the blueprint's structure intact and address its main gaps. A1 is mandatory; A2 raises the ceiling; A3 is an optional novelty upgrade.

## A1 — Remove outcome leakage from SRIS *(amends §3.1, node C)*
**What:** Redefine the clinical node `C` as **baseline-clinical only** — `age, grade, purity, Karnofsky`. Move `os_months` and `deceased/censored` out of the stalk; they become **endpoint-only** variables used in `survival_analysis.py`.
**Why:** As written, `C` contains survival, so `r_DC`/`r_RC` (and thus SRIS) are functions of survival — yet the headline claim is "SRIS predicts survival beyond covariates." That is circular and will sink the main result.
**How:** Drop the two outcome fields from `build_node_matrices`; verify SRIS no longer changes when survival is permuted.

## A2 — Promote external validation to a core aim *(amends §4.8, §5.5, §8)*
**What:** Move **CGGA replication** from "next 1 week, if feasible" to a primary result, and add it as an acceptance criterion: *"SRIS's prognostic effect replicates in CGGA."*
**Why:** Single-cohort prognostic scores are routinely rejected; cross-cohort replication is the strongest evidence a reviewer can get.
**How:** Harmonize CGGA to `clean_glioma_sheaf_table.csv`; **freeze** the discovery maps/weights on TCGA; recompute SRIS and re-run the Cox + group tests on CGGA. Use OT (already in §4.8) to align the two cohorts before transfer.

## A3 — Optional v2: patient-graph sheaf + prediction reliability *(new Phase 8, gated)*
**What:** A second sheaf — over the **patient similarity graph** (vertices = patients, edges = molecular kNN) — with learned O(d) restriction maps and **neural sheaf diffusion** for survival; emit each patient's local consistency radius as a **prediction-reliability flag** ("trust this prognosis or not").
**Why:** The 3-node modality sheaf is topologically thin; a reviewer may call SRIS "multi-view residuals in sheaf clothing." The patient-graph sheaf adds genuine ML novelty and a clinically actionable uncertainty output, differentiating the method.
**How:** Build only **after** v1 (the 3-node SRIS) clears its acceptance criteria. Keep stalk dim small, O(d)-constrained, nested CV — n≈420 is tight.

---

**One-line summary:** keep the blueprint as the execution plan; **do A1 now** (it's a correctness fix), **commit to A2** (it's the publishability lever), and **hold A3** as the stretch upgrade.
