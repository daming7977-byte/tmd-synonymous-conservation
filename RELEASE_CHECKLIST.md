# Final public-release sequence

1. Use the existing repository only: https://github.com/daming7977-byte/tmd-synonymous-conservation. Do not create a second repository. Before release, confirm the local working tree is clean and local `main` matches `origin/main`.
2. Use the existing Zenodo draft with reserved DOI `10.5281/zenodo.22705327`. Do not create a second Zenodo record and do not enable a separate GitHub-integrated archival record for this version.
3. When final release is authorized, change the existing GitHub repository from private to public. Then synchronize README, CITATION.cff, PRE_RELEASE_AUDIT.json and release metadata to the factual public state; regenerate repository hashes and commit the metadata-only release-state update. No scientific analysis rerun is authorized.
4. Create and push tag `v1.0.0` at the final reviewed release commit and publish the GitHub v1.0.0 release using `RELEASE_NOTES_v1.0.0.md`.
5. Generate the final archive from the reviewed `v1.0.0` tag with `git archive`, not Finder compression. Confirm the archive contains only tracked release files, excludes `.git`, backups, raw/third-party data and temporary files, and passes the internal `SHA256SUMS.txt` checks.
6. In the existing Zenodo draft for `10.5281/zenodo.22705327`, delete the draft-stage ZIP and upload the final archive generated from tag `v1.0.0`. Confirm title, creator, version, mixed/scoped licensing metadata, repository URL and file contents before publishing.
7. Publish that same Zenodo draft so `10.5281/zenodo.22705327` becomes registered. Confirm the DOI resolves and the public GitHub release is accessible.
8. Record the final public URL, tag/commit, DOI and archive hash. Update manuscript Data/Code Availability and the factual release-status wording, remove `PRE-RELEASE REVIEW`, regenerate the GBE submission PDF, and perform a final submission-only audit. Do not change scientific results or rerun analyses.

MIT applies only to author-owned code under `code/`; no blanket MIT license applies to the complete archive.
