"""
Lambda Tester for OT Cost Function
===================================
Tests different lambda values for:
  Cost(i,j) = ||s_i - s_j||^2 + lambda1 * SubtypePenalty(i,j) + lambda2 * |Age_i - Age_j|

Usage:
  python lambda_tester.py                          # uses defaults
  python lambda_tester.py --l1 0.4 --l2 0.02      # single run
  python lambda_tester.py --sweep                  # sweeps a grid of lambdas
  python lambda_tester.py --l1 0.3 --l2 0.05 --output results.txt
"""

import argparse
import math
import os
import sys
from itertools import product as iterproduct


# ── data loading ──────────────────────────────────────────────────────────────

def load_data(path="data.txt"):
    patients = []
    with open(path) as f:
        headers = f.readline().upper().split()
        for line in f:
            vals = line.split()
            if not vals:
                continue
            p = dict(zip(headers, vals))
            patients.append(p)
    return patients


# ── encoding (mirrors EncodedCases) ──────────────────────────────────────────

def encode_idh(p):
    """1 = mutant, 0 = wildtype, None = missing"""
    raw = p.get("IDH", "NA").strip().lower()
    if raw in ("na", "null", ""):
        return None
    return 1 if raw == "mutant" else 0


def encode_mgmt(p):
    """MGMT methylation: Methylated=1, Unmethylated=0, None=missing"""
    raw = p.get("MGMT", "NA").strip().lower()
    if raw in ("na", "null", ""):
        return None
    return 1 if raw == "methylated" else 0


def encode_egfr(p):
    """EGFR: amp=1, normal=0, None=missing"""
    raw = p.get("EGFR", "NA").strip().lower()
    if raw in ("na", "null", ""):
        return None
    return 1 if raw == "amp" else 0


def encode_months(p):
    """Overall survival in months, None if missing"""
    raw = p.get("MONTHS", "NA").strip()
    if raw in ("na", "null", ""):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def encode_grade(p):
    """Grade4: 1 if G4, 0 otherwise, None if missing"""
    raw = p.get("NEOPLASM_HISTOLOGIC_GRADE", p.get("NEOPLASMHISTOLOGICGRADE", "NA")).strip().upper()
    if raw in ("NA", "NULL", ""):
        return None
    return 1 if raw == "G4" else 0


def encode_age(p):
    """Age at diagnosis, None if missing"""
    raw = p.get("DIAGNOSISAGE", p.get("AGE", "NA")).strip()
    if raw in ("na", "null", ""):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def encode_subtype(p):
    """GBM vs LGG inferred from IDH status (wildtype -> GBM, mutant -> LGG)"""
    idh = encode_idh(p)
    if idh is None:
        return "unknown"
    return "GBM" if idh == 0 else "LGG"


# ── consistency (mirrors Consistency class) ───────────────────────────────────

def severity(grade4, months):
    if grade4 is None or months is None:
        return None
    return 1 if (grade4 == 1 or months <= 24) else 0


def r_DR(idh):
    if idh is None:
        return None
    return 1 - idh


def r_DC(idh, mgmt):
    if idh is None or mgmt is None:
        return None
    return 1 if (idh == 0 or mgmt == 0) else 0


def r_RC(egfr):
    return egfr


def e_DR(idh, egfr):
    r = r_DR(idh)
    if r is None or egfr is None:
        return None
    return abs(r - egfr)


def e_DC(idh, mgmt, grade4, months):
    r = r_DC(idh, mgmt)
    sev = severity(grade4, months)
    if r is None or sev is None:
        return None
    return abs(r - sev)


def e_RC(egfr, grade4, months):
    r = r_RC(egfr)
    sev = severity(grade4, months)
    if r is None or sev is None:
        return None
    return abs(r - sev)


def sheaf_vector(p):
    """Returns [eDR, eDC, eRC] for a patient dict, using None for missing."""
    idh   = encode_idh(p)
    mgmt  = encode_mgmt(p)
    egfr  = encode_egfr(p)
    months = encode_months(p)
    grade4 = encode_grade(p)

    edr = e_DR(idh, egfr)
    edc = e_DC(idh, mgmt, grade4, months)
    erc = e_RC(egfr, grade4, months)
    return [edr, edc, erc]


def sheaf_sq_diff(sv_i, sv_j):
    """||s_i - s_j||^2 — skips None dimensions."""
    total = 0.0
    counted = 0
    for a, b in zip(sv_i, sv_j):
        if a is not None and b is not None:
            total += (a - b) ** 2
            counted += 1
    return total if counted > 0 else None


# ── cost function ─────────────────────────────────────────────────────────────

def cost(pi, pj, l1, l2):
    """
    Cost(i,j) = ||s_i - s_j||^2
              + l1 * SubtypePenalty(i,j)
              + l2 * |Age_i - Age_j|
    Returns dict with total and each component.
    """
    sv_i = sheaf_vector(pi)
    sv_j = sheaf_vector(pj)

    sq_diff = sheaf_sq_diff(sv_i, sv_j)
    if sq_diff is None:
        sq_diff = 0.0

    subtype_i = encode_subtype(pi)
    subtype_j = encode_subtype(pj)
    subtype_pen = l1 if (subtype_i != subtype_j and
                         "unknown" not in (subtype_i, subtype_j)) else 0.0

    age_i = encode_age(pi)
    age_j = encode_age(pj)
    age_pen = l2 * abs(age_i - age_j) if (age_i is not None and age_j is not None) else 0.0

    return {
        "total":       sq_diff + subtype_pen + age_pen,
        "sq_diff":     sq_diff,
        "subtype_pen": subtype_pen,
        "age_pen":     age_pen,
    }


# ── output helpers ────────────────────────────────────────────────────────────

def patient_id(p):
    return p.get("PATIENTID", p.get("ID", "?"))


def print_cost_matrix(patients, l1, l2, out):
    n = len(patients)
    ids = [patient_id(p) for p in patients]
    col_w = 10

    header = f"{'':20}" + "".join(f"{pid[-col_w:]:>{col_w}}" for pid in ids)
    out.write(header + "\n")
    out.write("-" * len(header) + "\n")

    all_costs = []
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(None)
            else:
                c = cost(patients[i], patients[j], l1, l2)
                row.append(c["total"])
                all_costs.append(c["total"])
        matrix.append(row)

    # find min/max for annotation
    mn = min(all_costs) if all_costs else 0
    mx = max(all_costs) if all_costs else 1

    for i, row in enumerate(matrix):
        line = f"{ids[i][-20:]:20}"
        for j, val in enumerate(row):
            if val is None:
                line += f"{'  —  ':>{col_w}}"
            else:
                line += f"{val:>{col_w}.4f}"
        out.write(line + "\n")

    out.write("\n")
    out.write(f"  Range: {mn:.4f} (low) → {mx:.4f} (high)\n")
    return matrix


def print_best_matches(patients, l1, l2, out):
    n = len(patients)
    out.write(f"\n{'Patient':<22} {'Best match':<22} {'Cost':>8}  "
              f"{'||Δs||²':>8}  {'SubtypePen':>10}  {'AgePen':>8}\n")
    out.write("-" * 82 + "\n")

    for i in range(n):
        best_j, best_c = None, None
        for j in range(n):
            if i == j:
                continue
            c = cost(patients[i], patients[j], l1, l2)
            if best_c is None or c["total"] < best_c["total"]:
                best_j, best_c = j, c
        if best_j is None:
            continue
        out.write(
            f"{patient_id(patients[i])[-22:]:<22} "
            f"{patient_id(patients[best_j])[-22:]:<22} "
            f"{best_c['total']:>8.4f}  "
            f"{best_c['sq_diff']:>8.4f}  "
            f"{best_c['subtype_pen']:>10.4f}  "
            f"{best_c['age_pen']:>8.4f}\n"
        )


def run_single(patients, l1, l2, out, verbose=True):
    title = f"λ₁ = {l1:.4f}   λ₂ = {l2:.4f}"
    sep = "=" * max(len(title) + 4, 82)
    out.write(f"\n{sep}\n  {title}\n{sep}\n\n")

    if verbose:
        out.write("COST MATRIX\n\n")
        print_cost_matrix(patients, l1, l2, out)
        out.write("\nBEST MATCHES (lowest cost per patient)\n")

    print_best_matches(patients, l1, l2, out)


# ── sweep mode ────────────────────────────────────────────────────────────────

def sweep(patients, l1_values, l2_values, out):
    out.write("\nSWEEP SUMMARY — average cost of best matches\n")
    out.write(f"{'λ₁':>8}  {'λ₂':>8}  {'avg best cost':>14}  {'min best':>10}  {'max best':>10}\n")
    out.write("-" * 56 + "\n")

    results = []
    n = len(patients)
    for l1, l2 in iterproduct(l1_values, l2_values):
        best_costs = []
        for i in range(n):
            best = min(
                cost(patients[i], patients[j], l1, l2)["total"]
                for j in range(n) if j != i
            )
            best_costs.append(best)
        avg = sum(best_costs) / len(best_costs)
        mn  = min(best_costs)
        mx  = max(best_costs)
        results.append((l1, l2, avg, mn, mx))
        out.write(f"{l1:>8.4f}  {l2:>8.4f}  {avg:>14.4f}  {mn:>10.4f}  {mx:>10.4f}\n")

    # recommend the combo with the lowest average best-match cost
    best = min(results, key=lambda r: r[2])
    out.write(f"\n  Lowest avg best-match cost: λ₁={best[0]:.4f}, λ₂={best[1]:.4f}  "
              f"(avg={best[2]:.4f})\n")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_float_list(s):
    return [float(x) for x in s.split(",")]


def main():
    parser = argparse.ArgumentParser(description="Test lambda values for the OT cost function.")
    parser.add_argument("--data",   default="data.txt",  help="Path to data file")
    parser.add_argument("--l1",     type=float, default=0.3,  help="λ₁ subtype penalty")
    parser.add_argument("--l2",     type=float, default=0.05, help="λ₂ age penalty")
    parser.add_argument("--output", default=None,             help="Write results to this file")
    parser.add_argument("--sweep",  action="store_true",      help="Sweep a grid of lambda values")
    parser.add_argument("--l1-range", default="0.0,0.1,0.2,0.3,0.5,1.0",
                        help="Comma-separated λ₁ values for sweep")
    parser.add_argument("--l2-range", default="0.0,0.01,0.02,0.05,0.1,0.2",
                        help="Comma-separated λ₂ values for sweep")
    parser.add_argument("--limit",  type=int, default=None,
                        help="Only use first N patients (handy for quick tests)")
    args = parser.parse_args()

    # resolve data path relative to script location
    data_path = args.data
    if not os.path.exists(data_path):
        alt = os.path.join(os.path.dirname(__file__), args.data)
        if os.path.exists(alt):
            data_path = alt
        else:
            print(f"Error: cannot find data file '{args.data}'", file=sys.stderr)
            sys.exit(1)

    patients = load_data(data_path)
    if args.limit:
        patients = patients[: args.limit]

    out = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout

    try:
        out.write(f"Loaded {len(patients)} patients from '{data_path}'\n")

        if args.sweep:
            l1_vals = parse_float_list(args.l1_range)
            l2_vals = parse_float_list(args.l2_range)
            out.write(f"Sweeping λ₁ ∈ {l1_vals}\n")
            out.write(f"         λ₂ ∈ {l2_vals}\n")
            sweep(patients, l1_vals, l2_vals, out)
            out.write("\n--- Single run with selected defaults ---\n")
            run_single(patients, args.l1, args.l2, out, verbose=True)
        else:
            run_single(patients, args.l1, args.l2, out, verbose=True)

    finally:
        if args.output:
            out.close()
            print(f"Results written to '{args.output}'")


if __name__ == "__main__":
    main()
