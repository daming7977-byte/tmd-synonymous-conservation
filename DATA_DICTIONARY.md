# Data dictionary

TSV files retain archived numeric strings, including source precision. Blank cells mean not supplied/not applicable and must not be converted to zero. `source_file` / `source_location` provide provenance. `evidence_class`, `timing` and `interpretation_note` are editorial annotations, not new scientific results.

ST1: `N` is TMD count except gene companions (gene count); `clusters` or `genes` indicates gene clusters where supplied. `beta`, `se`, `ci_low`, `ci_high`, `p` are archived coefficient, standard error, 95% CI endpoints and P. `holm_p` applies only to the named three-factor family; `spearman_rho/p` are saved gene-companion statistics. TMD score scaling uses 5538 bridge TMDs; gene score scaling uses 1400 bridge genes. Shared temporal rows preserve audit precision and leave unreported SE blank.

ST2: `model`, `location`, `beta_per_1SD_evolution_score`, `cluster_robust_SE`, `CI_2.5`, `CI_97.5`, `p`, `Holm_p` are original fields. Spatial standardization uses 4505 TMDs and a distinct outcome; do not pool with ST1. Joint table `N_stacked` is not unique TMD count. Within-protein rows carry a discrepancy note.

S1 summaries: matched_N, Pearson_r, Spearman_rho, OLS_slope_reconstruction_on_author and OLS_intercept are saved validation results; residual_p05/p95 are residual percentiles, not CIs. S2: SMD_included_minus_excluded and multiplicity fractions are archived descriptive quantities. S3: expected_log2_shift is the profile-wide constant; beta_delta/se_delta/p_delta are saved numerical validation differences, not biological effects. Max direct codon rounding differences and canonical nonconstant residuals are distinct fields.
