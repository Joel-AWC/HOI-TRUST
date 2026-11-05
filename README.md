# HOI-TRUST

Explainable and uncertainty-aware machine learning for hospital-onset infection surveillance.

This release contains analysis source, recorded environment specifications and aggregate results
for the completed retrospective internal-validation run. It provides a research prototype;
external, genuinely temporal and prospective clinical validation have not been performed.
The compatibility partition is non-confirmatory because it informed earlier development.

## Contents

- `HOI_TRUST_pipeline.ipynb`: the original 52 code-cell sources, with saved outputs, cell metadata,
  execution counts, attachments and widget state removed.
- `results/`: aggregate numerical tables. Counts and metrics are conditional on the matched sample.
- `environment/`: the recorded direct dependency lock, full environment freeze and Python version.
- `scripts/regenerate_aggregate_figures.py`: plotting from the public aggregate CSVs only.
- `scripts/export_hoi_release.py`: the export script, when executed from its standalone `.py` file.
- `release_manifest.json`: run identity, source-export hashes, transformations and checks.
- `SHA256SUMS.txt`: checksums for the other release files.

The manifest's code-cell-source hash uses newline-joined cell sources. It is distinct from the
original notebook's embedded code-version hash, which uses the pipeline's own hashing convention.
The stripped notebook also has a different whole-file hash from the retained executed notebook.

## Source data and analysis execution

The source dataset is PhysioNet's *Predictors of Hospital Onset Infection*, version 1.0.0:
https://doi.org/10.13026/k70x-0m81. Credentialing, required training and the applicable data-use
agreement are required. Source records, individual predictions, models and private caches are
not included in this release.

To run the original analysis, use Python **3.13.16** and the recorded dependencies. Obtain authorised
source files and configure `HOI_TRUST_DATA_DIR` for a private project folder. Save a copy of the
notebook there under its expected name, `HOI_TRUST_submission_ready.ipynb`, and set
`HOI_TRUST_NOTEBOOK_PATH` to that copy before executing the configuration cell. Follow the run
instructions in the notebook, including the package-installation restart. Use
`FAST_MODE=False` and `HOI_TRUST_PROFILE=definitive` for the reported analysis. The notebook
checks source-release counts and data, code, environment and cache identities. Source-specific
fingerprints can differ on a different authorised release. Do not bypass those checks to mix runs.

The exporter does not require the current interpreter to match the modelling environment: it
copies or reconstructs the **recorded** environment and verifies its original hashes. It never
uses a new `pip freeze` as a substitute for the historical environment.

## Regenerate the public aggregate figures

From the extracted release root, with NumPy, pandas and Matplotlib available:

```bash
python scripts/regenerate_aggregate_figures.py
```

Seven PNG/PDF aggregate reproductions will be written to `figures/`. This command reads no
individual records, loads no fitted models and performs no model fitting or performance resampling.
These plotting scripts were added for the release; the original analysis notebook remains available
for the original plotting procedure. Layouts may differ from the submitted figure files.

| Manuscript item | Aggregate numerical sources in `results/` |
|---|---|
| Table 1 | `candidate_summary.csv`, `candidate_metrics_by_repeat.csv` |
| Table 2 and Figure 1 | `component_bootstrap_95ci.csv`, `paired_stack_vs_best_single_95ci.csv` |
| Table 3 | `probability_calibrator_apparent_fit.csv`, `v2_locked_operating_policies.csv` |
| Figure 2 | `strict_sensitivity_results.csv` (strict nested stack rows only) |
| Figure 3 | `capacity_curve_operating_patient.csv`, `capacity_curve_operating_environmental.csv`, `selective_curve_public_patient.csv`, `selective_curve_public_environmental.csv` |
| Figure 4 | `v2_legacy_uncertainty_metrics_NOT_CONFIRMATORY.csv` |
| Figure 5 | `legacy_within_pathogen_metrics_NOT_CONFIRMATORY.csv` |
| Figure 6 | `v2_locked_candidate_grouped_permutation_importance.csv` |
| Figure 7 | `dashboard_review_queue.csv`, `dashboard_uncertainty_composition.csv` |
| Main cohort/partition tables | `pathway_identity_summary.csv`, `analytical_partition_counts.csv` |
| Appendix calibration/subgroups/context sensitivity | `crossfitted_internal_calibration_by_repeat.csv`, `v2_legacy_subgroup_robustness_NOT_CONFIRMATORY.csv`, `pathogen_context_sensitivity.csv` |

The saved compatibility metrics CSV contains fresh 5%-capacity diagnostics; transferred locked
10%-capacity results are in `v2_legacy_capacity_metrics_NOT_CONFIRMATORY.csv`. Do not substitute
one for the other. Repeat-specific calibration estimates and candidate repeat means can differ
from the primary estimate formed by averaging repeated predictions by observation before resampling.

## Public aggregation and interpretation

The exporter selects named aggregate fields only. Small subgroup rows with fewer than 30 cases
or 30 controls are omitted. Dashboard data are exported as two margins: pathway/pathogen review
counts and pathway/uncertainty counts. Individual scores, score means and the original three-way
dashboard cells are omitted.

The selective-error curve is shared on an approximately 1%-spaced coverage grid with at least
30 retained observations and a minimum difference of 30 retained observations between adjacent
public points. Feasible locked-policy points are preserved. Per-point confidence thresholds are
not exported. Values are selected from the saved aggregate curve without interpolation; the
very small retained-fraction tail and fine-grained geometry of the submitted curve are omitted.
Exact fine-grained reconstruction requires authorised source access and the analysis notebook.
The 30-observation rule is an export convention, not a clinical performance guarantee or formal
proof of anonymisation.

Sensitivity estimates are single-repeat diagnostics without dedicated confidence intervals.
Permutation bars describe perturbation percentiles, not confidence intervals. Subgroups and
compatibility results are descriptive. Matched-sample probabilities are not population absolute risks.

The exporter verifies the saved run identifiers, the source notebook code, recorded environment
hashes, selected manuscript numerical anchors, CSV structure and public output checksums. This
does not independently establish the provenance of every CSV or reproduce the statistical analysis.
Inspect the release candidate before publication. Preserve the original restricted run separately.

## Version and citation

This package is prepared for release **v1.0.1**; the notebook's internal analysis version remains
3.1.0. After publishing and archiving the release, add its version-specific Zenodo DOI and correct
software-author metadata here and use that DOI in the manuscript. A release-candidate ZIP does
not establish that the files are already public.

## Licence

The author-generated code is MIT licensed, as in the existing repository. This does not replace
the source dataset's access conditions or grant permission to redistribute restricted records.
