# Account handoff and release sequence

1. Authenticate the author's existing GitHub account locally. Do not send passwords or tokens in chat. Confirm the actual destination owner and repository name in that session; no account identity is inferred from the manuscript author name.
2. Create a public repository from this directory only, not its enclosing private project. Do not upload parent audit folders, backup directories or legacy staging packages. A suitable proposed repository name is p5-evolutionary-functional-boundary; it is not a claimed existing URL.
3. After the real repository exists, add its actual URL to CITATION.cff and README, regenerate the current file hashes and archive, and commit the reviewed content. Do not assign a blanket MIT license to excluded/noncode content.
4. Sign into the existing Zenodo account and enable the new GitHub repository, if using automatic integration, before publishing the GitHub release. Account sign-in and authorization require the user. Alternatively upload the reviewed release archive manually. Use the approved creator Ming Li; leave ORCID absent unless verified.
5. Tag and publish v1.0.0 with RELEASE_NOTES_v1.0.0.md. Confirm the resulting Zenodo record and its version-specific DOI actually exist. Ensure archival license metadata states the mixed/scoped rights accurately and does not imply MIT for the full archive. No .zenodo.json with a misleading blanket license is supplied.
6. Record actual URL, tag/commit, version DOI and release-file hashes. Update manuscript Data/Code Availability and the locked Methods' release-status sentence only factually; regenerate the submission PDF and perform the final submission audit. No analysis rerun is authorized.

Public release is approved, but this local preparation does not bypass user login. Official integration instructions: https://help.zenodo.org/docs/github/ .
