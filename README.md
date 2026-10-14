# RSV Single-cell Resource

The frozen resource contains 2,474,742 cell UIDs across 15 source labels. The reviewed compilation retains every project annotation row, 407 comparison eligibility records and 45 focused infection/control candidates (42 CONDITIONAL, 3 HOLD). Candidate/Unknown identities and review reasons are retained. Detailed clinical values and exact ages from 22 original field names remain withheld; NULL may mean intentional privacy withholding. Original source terms and reference-overlap limitations remain applicable.

At the recorded closure check, 65/65 files and 14/15 source H5ADs covering 2,360,607 cells passed server availability, byte-size and MD5 verification. The two previously pending files were replaced within their original private entries and are now available with matching byte size and server MD5; there are no duplicate registered files. The five-file full compilation is verified and was actually downloaded through the owner account, rehashed and reopened. All data entries remain private/unpublished. Confidential reviewer access is managed separately from this public code repository; its current availability and download validation are recorded in the confidential submission access report. No project DOI has been created. A private collection may organize the approved entries without itself granting anonymous download access. Owner authentication success does not demonstrate anonymous access.

Sanger/Yoshida has 114,135 frozen cells in the compilation. The original authors provide an official portal H5AD download, and the official object was tested with unauthenticated HEAD and a bounded range read. Its current multipart ETag matches the frozen input; a different portal whole-file MD5 is recorded. This supports original source acquisition. Direct redistribution of the project H5AD and canonical source mask/map remains unconfirmed. The compilation retains non-expression frozen selection, gene namespaces and generic reconstruction maps; it contains no Sanger expression counts. See config/sanger_reconstruction.json for the conditional third-party reconstruction dependencies and its explicit NOT_RUN result. Earlier internal reconstruction retained 114,135 cells across 54 capture groups with zero count or measurement-mask differences; this is not independent external validation.

After obtaining the reviewed files through an approved data access route, assemble `release/expression/`, `release/features/` and `release/curated_metadata/`. Extract the compilation tar into the last directory. The currently available code is read-only by default; source acquisition requires an explicit --download flag.

```bash
conda env create -f environment.yml
conda activate rsv-publication-reader
python scripts/download_sources.py --config config/verified_source_downloads.tsv --output-root ./original_sources
python scripts/load_resource.py --release-root ./release --query comparisons --focused --limit 50
python scripts/load_resource.py --release-root ./release --query identity --source CD1c_DC --limit 3
python scripts/load_resource.py --release-root ./release --query cells --source Berg_ALI --limit 2
python scripts/validate_resource.py --release-root ./release --output ./validation
```

The full source table covers 15 labels and 57 original acquisition objects. The list-only downloader was tested without downloading those full objects. Reconstruction CLI help was tested; expression reconstruction/export was NOT_RUN during closure. The validation command needs the three fixed source H5ADs, masks and `examples/FROZEN_MARKER_EXPECTATIONS.tsv` from the approved release; these data are not GitHub code assets. A missing expression file is reported separately from an available metadata/annotation row.

See DATA_ACCESS.md, DATA_DICTIONARY.tsv and config/resource_delivery.json for versions and access scope. Software/documentation use the MIT License. This does not relicense third-party data. Candidate labels and source/reference consistency do not constitute independent accuracy validation. Author facts and formal resource citation remain pending.
