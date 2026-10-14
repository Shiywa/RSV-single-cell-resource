# RSV single-cell resource

Source-specific expression, unified cross-study metadata, versioned annotation, and provenance-aware access. This repository contains reading, verification, acquisition and reconstruction code plus source configuration and usage examples. The software is released under the MIT License; see LICENSE_SCOPE.md for the distinction between software and data permissions. Author and publication facts remain pending in CITATION_PENDING.md.

## Data availability

The complete internal collection contains 2,474,742 cells across 15 source labels. The source-minimized deposit scope contains 818,148 cells from Berg_ALI, Qin_cotton_rat and HLCA in three source-specific H5ADs. Exact integer counts, original gene namespaces and source-level measurement masks are preserved. These scopes have different access and permission conditions.

Figshare authentication and private upload testing have succeeded. Permitted source-minimized assets are being archived as unpublished private drafts. No public project data download URL or project DOI is available. The repository does not contain expression matrices, clinical metadata, credentials or the complete project annotations. Original publicly available source-file URLs and checksums are listed in config/verified_source_downloads.tsv.

## Installation and use

Create the environment specified in environment.yml with conda or micromamba, then use the commands below with a legitimately obtained release directory.

```bash
conda env create -f environment.yml
conda activate rsv-publication-reader
python scripts/load_resource.py --release-root /path/to/authorised/release --query cells --source HLCA --limit 10
python scripts/validate_resource.py --release-root /path/to/authorised/release --output /path/to/validation-output
python scripts/download_sources.py --config config/verified_source_downloads.tsv --output-root /path/to/downloads
```

The final command only lists the original source downloads. Add --download explicitly to acquire those sources under their applicable terms. Reconstruction is a separate operation requiring authorized frozen selections and feature maps; it is not run by installing or reading this repository.

With an authorized complete release, additional queries include `--query comparisons --focused`, `--query identity --source CD1c_DC`, and `--query trace --cell <CellUID>`. Complete-project queries return BLOCKED_BY_RELEASE_SCOPE in the source-minimized candidate. `validate_resource.py` checks actual H5AD/count/mask/metadata joins; optional forbidden roots enforce local isolation. `reconstruct_h5ad.py` wraps the verified sparse native-to-source export. The reused `rebuild_original_sources.py` verifies legitimate source files against authorized frozen selection/maps/native objects; its historical Berg skip is disclosed and the existing full Berg reconstruction receipt remains part of the separately controlled provenance. Missing restricted dependencies cannot be reconstructed lawfully from a repository URL alone.

`download_sources.py` lists confirmed source-file URLs by default; downloading requires an explicit flag and source-compatible terms. No acquisition ran during publication packaging. Other12 sources require asset-specific rights and properly delivered selection/feature mapping/annotation dependencies. Fields in source_manifest describe15 study labels, with one HLCA primary object rather than11 view copies. Independent feature namespaces and exact int32 counts are retained.
