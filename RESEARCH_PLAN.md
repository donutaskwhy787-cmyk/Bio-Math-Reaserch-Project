# PINGO — Research Protocol & Execution Checklist

**Working title:** *A sheaf-theoretic patient-similarity framework localizes molecular–outcome discordance and prognostic unreliability in diffuse glioma*

**One-sentence thesis:** Using routine molecular markers we build a patient-similarity graph, equip it with a **cellular sheaf**, and show that the survival signal's **obstruction to gluing** (sheaf cohomology / Laplacian energy) defines *molecular–outcome discordance* and, crucially, *localizes where prognostic prediction is fundamentally unreliable* — a reproducible signal discovered in TCGA and validated in CGGA, with prediction by **neural sheaf diffusion**.

**Two mathematical pillars (PINGO):** **(1) Optimal transport** aligns the two cohorts' representation spaces so the model transfers; **(2) cellular sheaf theory** defines discordance, prediction, and per-patient uncertainty. Decided scope: **Tier B** — sheaf cohomology *defines discordance* AND neural sheaf diffusion *does prediction with uncertainty*. Team confirmed strong on both topology/algebra and GNNs.

**Foundational references to follow exactly (don't improvise the formalism):** Hansen & Ghrist, *Toward a Spectral Theory of Cellular Sheaves* (2019); Curry, *Sheaves, Cosheaves and Applications* (thesis); Robinson, *Topological Signal Processing* (consistency radius / data fusion); Bodnar, Di Giovanni, Chamberlain, Liò, Bronstein, *Neural Sheaf Diffusion* (NeurIPS 2022) + its open-source code.

**Cohorts (decided):**
- **Discovery:** TCGA pan-glioma (Ceccarelli et al., *Cell* 2016) — your `data.txt`, 420 patients.
- **Validation:** CGGA `mRNAseq_693` (primary) + `mRNAseq_325` (second replication).

**Why this can publish:** the dataset is heavily mined, so we do *not* claim new subtypes. The contributions are (1) a reproducible mixed-type similarity + OT-alignment **method**, and (2) a cross-cohort-replicated **catalog of discordant patients** with their distinguishing biology.

---

## Core design principles (read first — everything below serves these)

1. **Non-circularity.** Survival (`Months`, `Status`) is NEVER an input to the distance. It is the thing we explain/predict. Violating this invalidates the whole paper.
2. **Two-tier features.** A *transferable* panel (present in both cohorts) defines the distance that must replicate; a richer *discovery-only* panel only *describes* the discordant patients.
3. **Adjust before you're surprised.** "Outlier" means *unexpected given what we already know* (age, IDH, MGMT, 1p/19q, grade). Raw distance-outliers are not publishable; Cox-residual outliers are.
4. **Replication is mandatory.** A single-cohort outlier list will be rejected. CGGA replication is the spine of the paper.
5. **Everything is reproducible.** Seeded, scripted, version-controlled, one-command rerun.

---

## PHASE 0 — Project setup & reproducibility
- [ ] **Create a clean repo structure**: `data/raw/`, `data/processed/`, `src/`, `results/`, `figures/`, `manuscript/`. Keep raw data read-only.
- [ ] **Pin the environment**: Python ≥3.11; `numpy`, `pandas`, `scikit-learn`, `scipy`, `lifelines` (survival/Cox), `gower` or custom Gower, `POT` (OT), `matplotlib`/`seaborn`, `statsmodels`. Freeze to `requirements.txt`/`environment.yml`.
- [ ] **Set a global random seed** and record it; every stochastic step (kNN ties, bootstraps, imputation) reads it.
- [ ] **Make a single driver** (`run_all.py` or a Makefile) that regenerates every result/figure from raw inputs. Reviewers and you-in-3-months both need this.
- [ ] **Start a lab notebook / decisions log** capturing every analytic choice and its justification (feeds the Methods section).

## PHASE 1 — Data acquisition & harmonization
- [ ] **Lock the TCGA discovery table** from `data.txt` (already in hand). Snapshot it to `data/raw/` and never edit in place.
- [ ] **Download CGGA** `mRNAseq_693` and `mRNAseq_325`: the clinical/annotation TSV (IDH, 1p/19q, MGMT, grade, age, sex, OS, censoring) + the expression matrix (for RNA-derived features). Source: https://www.cgga.org.cn/download.jsp
- [ ] **Build a harmonization map** — a documented dictionary translating each cohort's coding to a shared schema. Examples to resolve:
  - IDH: TCGA `Mutant/WT` ↔ CGGA `Mutant/Wildtype`.
  - 1p/19q: TCGA `IDH/codelsubtype` (codel/non-codel) ↔ CGGA `Codel/Non-codel`.
  - MGMT: TCGA `Methylated/Unmethylated` ↔ CGGA methylation status (confirm same assay basis).
  - Grade: standardize to integer 2/3/4 (handle WHO-2021 vs older labels).
  - Survival: TCGA `Months` ↔ CGGA OS (confirm units = months; convert days if needed); `Status` ↔ CGGA censoring (1=event/dead, 0=censored/alive).
- [ ] **Define the transferable feature panel** (must exist in both): `IDH`, `MGMT`, `1p19q-codel`, `grade`, `age`, and RNA-derived `EGFR-expr`, `ESTIMATE-immune`, `ESTIMATE-stromal`, `TERT-expr` (computed identically per cohort — see Phase 3).
- [ ] **Define the discovery-only panel** (TCGA, characterization only): `Chr7gain/Chr10loss`, `Chr19/20co-gain`, `ATRXstatus`, `BRAFV600Estatus`, `TMB`, `Percentaneuploidy`.

## PHASE 2 — Cleaning, QC & missing data
- [ ] **Type-cast on ingest.** Current code loads everything as strings then does numeric comparisons — fix at the loader so numbers are numbers, categoricals are categoricals.
- [ ] **Drop non-features:** identifiers (`PatientID`, `SampleID`), constants (`CancerType`, `CancerTypeDetailed`).
- [ ] **Drop >60%-missing columns:** both telomere-length fields (83%), `TelomereMaintenance` (69%), `TERTpromoterstatus` (69%).
- [ ] **Drop collinear duplicates:** `MutationCount` (keep `TMB`), `ESTIMATEcombinedscore` (keep immune+stromal), `TERTexpressionstatus` (keep `TERTexpression`).
- [ ] **Impute 20–45%-missing kept features** (ATRX, BRAFV600E, TMB, ESTIMATE, TERT-expr): median for continuous, mode for categorical, **plus a binary `was_imputed` flag per feature** so models know it wasn't measured. (Your code already TODOs this.)
- [ ] **Confound handling (decided):** `AbsolutePurity` is a *technical* covariate — use it to sanity-check/correct expression-derived features, not as a similarity feature; `Sex` enters only as a Cox covariate, not the distance; `KarnofskyPerformanceScore` excluded from distance (semi-outcome, 33% missing).
- [ ] **QC report:** per-cohort missingness table, value distributions, and a flag list of suspicious samples (extreme purity, contradictory IDH/subtype labels).

## PHASE 3 — Feature engineering & the cost function
- [ ] **Standardize continuous features** (z-score within each cohort) BEFORE combining. This is the fix for the dominance bug we measured: today age contributes ~304 and survival ~141 to a cost where flipping IDH moves it 0.5 — i.e., the molecular weights are currently ignored.
- [ ] **Compute RNA-derived features identically in both cohorts:** run the **ESTIMATE** algorithm on each expression matrix → immune/stromal scores; extract `EGFR` and `TERT` expression (log2, same normalization). This is what makes them comparable across TCGA and CGGA.
- [ ] **Encode types correctly:** grade as **ordinal** (2<3<4, so G2↔G4 costs more than G2↔G3 — the current 0/1 flip loses this); nominal markers (IDH, MGMT, codel, BRAF, ATRX) as categorical mismatch.
- [ ] **Implement the cost function as Gower distance** over the panel:
  - per continuous feature: scaled absolute difference in [0,1];
  - per ordinal: normalized rank difference;
  - per nominal: 0 if equal, 1 if different;
  - aggregate = weighted mean over features **both patients actually have** (pairwise-complete), normalized by the contributing count.
- [ ] **Fix the missing-value branch** in the old `omicsDistance` (its NA condition is logically unreachable) — Gower's pairwise-complete averaging replaces it cleanly.
- [ ] **Set weights explicitly and justify them** (e.g., IDH/codel highest as WHO-defining). Default to equal weights for the main analysis; treat the chosen weights as a hypothesis to stress-test in Phase 7.
- [ ] **Output:** an N×N cost/similarity matrix per cohort (no OT yet — for discordance we use the matrix directly).

## PHASE 4 — Similarity graph, OT alignment & the cellular sheaf
- [ ] **Build the within-cohort distance matrix** from the cost function (one table in → full pairwise matrix out; diagonal 0).
- [ ] **Build the patient graph G=(V,E):** vertices = patients, edges = molecular k-NN (sweep k ∈ {5,10,15,20}; pick by validation stability, not by the result you want). This graph is the **base space** of the sheaf.
- [ ] **OT's job (decided) — cross-cohort alignment:** entropic OT / Sinkhorn (`POT`) aligns the TCGA and CGGA representation spaces (domain adaptation) so the sheaf's restriction maps and obstruction modes transfer across cohorts. OT also gives a principled metric to compare the two cohorts' sheaf spectra.
- [ ] **Equip G with a cellular sheaf F (the core construct):**
  - **Stalks:** vertex stalk `F(v)=ℝ^d` (small d, e.g. 2–8), a latent representation of patient *v* initialized from their **molecular** feature vector; edge stalk `F(e)=ℝ^d`.
  - **Restriction maps:** for edge `e=(u,v)`, `F_{v◁e}: F(v)→F(e)`, learned/derived from molecular features. **Decided:** constrain to **O(d) orthogonal** maps (norm-preserving, fewer parameters, identifiable at small n, interpretable) per the O(d)-sheaf variant.
  - **Coboundary & Laplacian:** `(δx)_e = F_{u◁e}x_u − F_{v◁e}x_v`; **sheaf Laplacian** `L_F = δᵀδ` (block diag `L_vv = Σ_{e∋v} F_{v◁e}ᵀF_{v◁e}`, off-diag `L_uv = −F_{u◁e}ᵀF_{v◁e}`). The trivial sheaf (d=1, identity maps) recovers the ordinary graph Laplacian — that's the baseline to beat.
- [ ] **GO/NO-GO diagnostic (do this early):** confirm the cohort carries **non-trivial sheaf structure** — i.e., the survival signal has meaningful high-energy (obstruction) components and `H⁰` is not degenerate. If everything glues trivially, there is no sheaf story → fall back to the Phase-5 Cox framing. Decide this gate before investing in the full model.

## PHASE 5 — Discordance scoring (the core analysis)
**Three definitions, reported together for robustness — the sheaf one is the novelty, the others are baselines it must beat.**
- [ ] **(Baseline 1) Cox model-based.** Fit Cox PH on known prognostic factors only (age, IDH, MGMT, 1p/19q, grade) — "what we already know." Check PH assumptions (Schoenfeld). Discordance = martingale/deviance residual.
- [ ] **(Baseline 2) Neighborhood-based.** Difference between a patient's outcome and their molecular k-NN neighbors' outcomes (observed vs. neighbor-average survival / neighbor KM).
- [ ] **(Novel) Sheaf-cohomological discordance.** Build the sheaf restriction maps from **molecular features only** (geometry is molecular → keeps survival out, non-circular). Treat the **survival/risk signal as a 0-cochain** and decompose it in the **sheaf-Laplacian spectrum**:
  - the low-energy / `H⁰` component = the part of survival variation that *glues* smoothly with molecular neighborhoods (expected);
  - the **high-energy / obstruction component = the part that refuses to glue = discordance.**
  - Per-patient score = that patient's local contribution to the Dirichlet energy `xᵀL_F x` / their **consistency radius**. This is the sheaf analog of graph-signal non-smoothness, and it is the paper's central definition.
- [ ] **Define discordant patients** by a pre-registered threshold (e.g., top & bottom deciles), set BEFORE looking at biology. Separate the two clinically meaningful tails:
  - *Favorable-discordant:* poor-prognosis molecular profile, unexpectedly long survival ("exceptional survivors").
  - *Unfavorable-discordant:* good-prognosis profile, unexpectedly short survival.
- [ ] **Agreement check:** quantify overlap of the three discordance rankings (Cox vs neighborhood vs sheaf). The sheaf must add information, not merely reproduce the Cox residual (else drop it — see Phase 7).
- [ ] **Stability check:** bootstrap the whole pipeline; report how often each patient is flagged (only stable flags enter the catalog).

## PHASE 5b — Prediction by neural sheaf diffusion (the prediction track)
- [ ] **Model:** discretized sheaf diffusion `ẋ = −L_F(t)x` with **learned O(d) restriction maps** (Bodnar et al.; adapt their open-source code). Final stalk representations → a **survival head** (Cox partial-likelihood / DeepSurv-style, or discrete-time hazard).
- [ ] **Built-in uncertainty (the oncology payoff):** output, per patient, both a risk prediction **and** the local **consistency radius / residual sheaf energy** as a *prediction-reliability* flag — "trust this prognosis or not." This is the clinically actionable contribution; build it as a first-class output, not an afterthought.
- [ ] **Small-n safeguards:** low stalk dimension d, O(d) constraint, dropout/weight-decay, early stopping, **nested cross-validation**; never tune on CGGA.
- [ ] **Restriction-map interpretation:** inspect/constrain learned maps for biological/structural meaning (orthogonal rotations between molecular contexts) — required so the model isn't "just a fancy GNN."
- [ ] **Outputs:** risk score, discordance/obstruction score, and reliability flag — all per patient, in both cohorts.

## PHASE 6 — Cross-cohort validation
- [ ] **Freeze the discovery model** (features, weights, k, thresholds) on TCGA. No peeking at CGGA while tuning.
- [ ] **Apply the frozen rule to CGGA** (after OT alignment) and recompute discordance.
- [ ] **Replication endpoints:**
  - Does the discordance score stratify survival in CGGA (log-rank, KM curves)?
  - Do the *characteristics* of discordant patients (Phase 8) recur in CGGA?
  - Quantify transfer (e.g., C-index of the discordance score in CGGA vs. discovery).
- [ ] **If it fails to replicate:** report honestly and pivot to a methods-focused paper (the framework + negative result is still publishable if rigorous). Decide this rule *now*, not after seeing the data.

## PHASE 7 — Statistical rigor & baselines
- [ ] **Baselines to beat** (reviewers WILL ask "why not the simple thing?"):
  - WHO molecular subtype alone as the prognostic stratifier.
  - Plain standardized-Euclidean distance instead of Gower.
  - Cox model with no neighborhood term.
  Show your framework adds discriminatory value (ΔC-index, likelihood-ratio test).
- [ ] **Sheaf-specific ablations (MAKE-OR-BREAK — the sheaf must earn its place):**
  - **Trivial sheaf vs. learned sheaf:** sheaf Laplacian vs. ordinary graph Laplacian; **Neural Sheaf Diffusion vs. GCN/GAT** and a plain Cox/MLP — on both prediction (C-index) and discordance recovery.
  - **Restriction-map family:** O(d)-orthogonal vs. diagonal vs. general; justify the choice empirically.
  - **Stalk dimension d sweep**; show non-triviality (d>1 helps) without overfitting.
  - **Sheaf vs. Cox discordance:** the sheaf score must capture discordant patients the Cox residual misses (incremental value), or be dropped honestly.
- [ ] **Multiple-testing control** (Benjamini–Hochberg) anywhere you scan many features/patients.
- [ ] **Sensitivity analyses:** vary feature weights ±50%, k, distance metric (Gower vs OT vs Euclidean), and imputation strategy. Conclusions must survive; report where they don't.
- [ ] **Power/sample-size reality:** outliers are rare — state effective n for each tail and avoid over-claiming on tiny subgroups.
- [ ] **Pre-register** the analysis plan (even informally, internally dated) to defend against "you tried 100 things."

## PHASE 8 — Biological characterization & interpretation
- [ ] **Profile the discordant tails using the discovery-only panel** (Chr7/10, Chr19/20, ATRX, BRAF, TMB, aneuploidy) + the methylation/RNA subtype labels you held out as validation: what distinguishes exceptional survivors from their molecular neighbors?
- [ ] **Test enrichment** of any marker/subtype in each tail vs. the cohort (Fisher/χ², BH-corrected).
- [ ] **Sanity-check artifacts:** confirm discordant ≠ "lots of imputed features" and ≠ purity/assay extremes.
- [ ] **Clinical narrative:** translate findings into a hypothesis (e.g., "a subset of IDH-wildtype patients with X carry unexpectedly favorable outcomes") — the part clinicians and reviewers remember.

## PHASE 9 — Figures & tables
- [ ] **Fig 1** — study/workflow schematic (cohorts, feature tiers, method).
- [ ] **Fig 2** — the similarity space (UMAP/MDS of the distance matrix), colored by subtype, with discordant patients highlighted.
- [ ] **Fig 3** — discordance score vs. survival; KM curves for tails in TCGA **and** CGGA (the money figure: replication).
- [ ] **Fig 4** — biological characterization of the discordant tails.
- [ ] **Table 1** — cohort characteristics (TCGA vs CGGA).
- [ ] **Table 2** — the discordant-patient catalog (de-identified) with scores and key markers.
- [ ] **Supplement** — missingness, sensitivity analyses, baseline comparisons, full feature dictionary.

## PHASE 10 — Manuscript & submission
- [ ] **Methods first** (it's mostly written if the decisions log is kept) — distance, OT alignment, Cox, thresholds, validation, seeds.
- [ ] **Reporting standards:** follow REMARK (prognostic-marker studies) / TRIPOD (prediction models) checklists — reviewers expect them.
- [ ] **Data & code availability:** public repo + DOI (Zenodo). TCGA/CGGA are public, so this is easy and expected.
- [ ] **Author contributions & the shared planning doc** — reconcile with the group (the `CsPlanningDoc.txt` link).
- [ ] **Target journals (tiered):** methods/応用 venues like *Bioinformatics* / *GigaScience* / *BMC Bioinformatics* (method emphasis), or *Neuro-Oncology* / *Acta Neuropathologica Communications* (clinical emphasis) if replication is strong.

---

## Code refactor checklist (current `Main.py` issues to clear)
- [ ] **Remove `ot.solve_sample(c, c)`** — it re-derives a distance from the cost matrix and self-matches; it's the original bug.
- [ ] **Stop loading all values as strings** — type-cast at ingest (breaks `>= 0.3`, `< 24` comparisons today).
- [ ] **Fix `omicsDistance` NA logic** — replace with Gower pairwise-complete averaging.
- [ ] **Grade handling** — replace string-slice `[1:2]` parsing and 0/1 penalty with proper ordinal encoding.
- [ ] **De-hardcode** filenames/columns in `data`, `OTOutputs`, and the `rain…` dump function; parameterize the feature set.
- [ ] **Standardize features** before combining (the dominance fix).
- [ ] **Delete/retire** unused scaffolding (`Consistancy`, half-built `EncodedCases` methods) or finish them deliberately.
- [ ] **Add tests** on the distance (symmetry, zero diagonal, known toy cases) and a seeded end-to-end smoke run.

## Risk register (decide mitigations up front)
- **No replication in CGGA** → pre-committed pivot to methods paper (Phase 6).
- **Missing data drives "outliers"** → imputation flags + artifact checks (Phases 2, 8).
- **Cross-cohort feature mismatch** → two-tier design + OT alignment (Phases 1, 4).
- **Over-fitting/over-claiming on rare tails** → bootstrap stability, power statement, multiple-testing control (Phases 5, 7).
- **"It's just GBM/IDH-wt"** → Cox adjustment + baseline comparisons (Phases 5, 7).
- **Sheaf is decoration (adds nothing over a GNN/Cox)** → mandatory ablations (Phase 7); pre-committed rule to drop it if the trivial sheaf ties.
- **No non-trivial sheaf structure in the data** → early GO/NO-GO diagnostic (Phase 4) before investing in Tier B.
- **Restriction maps unidentifiable at small n** → O(d) constraint, low d, regularization, nested CV (Phases 4, 5b).
- **Math reviewers find formalism errors** → follow Hansen–Ghrist / Bodnar exactly; have the math member own the appendix proofs.

## Critical path (do in this order)
1. CGGA download + harmonization map (Phase 1) — unblocks everything cross-cohort.
2. Clean loader + Gower distance + standardization (Phases 2–3) — the corrected engine.
3. Patient graph + **sheaf GO/NO-GO diagnostic** (Phase 4) — decides whether Tier B is viable *before* heavy investment.
4. Cox/neighborhood/**sheaf** discordance on TCGA (Phase 5) — the core result.
5. **Neural sheaf diffusion** predictor + uncertainty (Phase 5b).
6. OT alignment + CGGA replication (Phases 4, 6) — the make-or-break.
7. Baselines + **sheaf ablations**, sensitivity, biology, figures (Phases 7–9).
8. Write-up (Phase 10).

---

## APPENDIX — Sheaf framework, formal summary (for the Methods/Supplement)
**Base space.** Patient k-NN graph `G=(V,E)` from the Gower distance.

**Cellular sheaf `F`.** Stalks `F(v)=F(e)=ℝ^d`. For `e=(u,v)`, orthogonal restriction maps `F_{u◁e}, F_{v◁e} ∈ O(d)` (learned, or derived from molecular features).

**Coboundary / Laplacian.** `δ: C⁰=⊕_v F(v) → C¹=⊕_e F(e)`, `(δx)_e = F_{u◁e}x_u − F_{v◁e}x_v`. Sheaf Laplacian `L_F = δᵀδ ⪰ 0`.

**Cohomology & energy.** `H⁰(G;F)=ker L_F` = global sections (assignments consistent across every edge). Dirichlet energy `E(x)=xᵀL_F x = Σ_{e=(u,v)} ‖F_{u◁e}x_u − F_{v◁e}x_v‖²`.

**Discordance (decided definition).** With molecular-only restriction maps, lift the survival/risk signal to a cochain `s`; decompose `s` over the eigenbasis of `L_F`. Low-energy mass = neighborhood-consistent (expected) survival; **high-energy / obstruction mass = discordance.** Per-patient discordance = local term `Σ_{e∋v} ‖F_{u◁e}s_u − F_{v◁e}s_v‖²` (≈ Robinson consistency radius at `v`).

**Prediction (decided model).** Neural sheaf diffusion `x(t+1)=x(t) − σ(L_{F(t)} (I⊗W₁) x(t))W₂`, learned `O(d)` maps, survival head with Cox partial-likelihood loss; per-patient consistency radius emitted as a **prediction-reliability** score.

**Cross-cohort transfer.** OT aligns TCGA and CGGA stalk spaces so `F`'s maps/spectrum are comparable; report transfer of both prediction (C-index) and obstruction patterns.

**Non-triviality checks.** dim `H⁰`, energy spectrum, and the Phase-7 ablations vs. the trivial (graph-Laplacian) sheaf establish that the sheaf structure is real and necessary.
