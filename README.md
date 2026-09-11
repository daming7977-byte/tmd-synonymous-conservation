# Evolutionary-functional membrane-protein boundary study

Reproducibility code and archived summary results for the combined human membrane-protein study. Author: Ming Li. Scientific analysis is closed. This repository candidate is prepared for version v1.0.0; no GitHub URL, published release or DOI is claimed until actually issued.

## Contents

- code/: project-authored scientific/validation scripts and the figure renderer; MIT applies only here.
- results/: ST1 with metadata-only evidence/timing alignment to the frozen submission table, plus unchanged ST2 and joint-test tables; numerical/statistical values are unchanged.
- source_summaries/: unchanged small archived descriptive/technical summaries, not raw or third-party author tables.
- manifests/: historical local contracts and selected freeze digests, with private path prefixes removed where needed.
- PROVENANCE.tsv and PATH_EDIT_LOG.tsv: original file digests and release-copy differences.
- ENVIRONMENT.md, DATA_DICTIONARY.md and SCRIPT_TO_RESULT_MAP.tsv: definitions, recorded environments and script/output roles.

## Reproduction limits and safe use

This is an archival code package, not a tested one-command portable workflow. Required external data/reference trees and some historical dependencies are not bundled. Scripts can write beside themselves or into an external project tree. Do not run them in a frozen archive. No scientific script or model was executed for this release.

Path-only release changes replace historical home-directory discovery and workstation tool paths with explicitly supplied environment values. Set P5_DATA_PARENT to the parent of the historical external data-tree folders; P5_PROJECT_ROOT to the human MPT project root; P5_FASTQ_DUMP and P5_MAP_ENV to the relevant installed tools; P5_MPLCONFIGDIR to a writable plotting cache and P5_RENDER_OUT to a separate rendering output directory. The remaining relative directory structure must match the source contracts; this package does not silently rebuild missing references. These names have no populated private defaults. Validation scripts retain their original relative-location assumptions and must be staged into the documented historical layout before any separately authorized execution.

Scientific formulas, coefficients, thresholds, windows, populations and analysis timing have not been changed. Original source hashes are distinguished from sanitized-release hashes. Syntax-only checks were performed; no inference/reconstruction execution is claimed. Historical environment gaps remain documented.

## Evidence boundaries

Primary MPT support was not established; adjusted evidence remains secondary. Profiles are not independent replications. Temporal localization and machinery specificity were not established. Historical spatial within-protein findings remain supporting with their specification–implementation discrepancy. Normalization restoration shifts the intercept without changing primary score inference. See tables and notes rather than interpreting isolated P values.

## Data and license

Public input datasets: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297497 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300977 . Obtain excluded inputs from their providers. See REDISTRIBUTION_EXCLUSIONS.md and LICENSE_SCOPE.md. MIT approval covers author-owned code only and is not a grant over third-party material.

## Release and archival status

CITATION.cff records the private GitHub repository URL and the reserved Zenodo DOI 10.5281/zenodo.22705327. No public GitHub release or registered DOI is claimed yet. The DOI will be registered when the Zenodo record is published. Follow RELEASE_CHECKLIST.md. No separately closed extension results are included.

Full manuscript freeze identifiers are retained in manifests/MANUSCRIPT_FREEZE_HASHES.md. This editorial release preparation did not execute scientific code or alter archived results.
