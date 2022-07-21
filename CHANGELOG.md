# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Changed
- Fastqingress metadata map

## [v0.3.2]
### Changed
- Improved `--help` command
- Improved robustness of `fastq_ingress` utility
- Rebuild docker container to include `longphase` v1.2
### Fixed
- Set `out_dir` option format to ensure output is written to correct directory on Windows

## [v0.3.1]
### Fixed
- `mosdepth` step will not return a blank result if a BED is not provided to the workflow

## [v0.3.0]
### Changed
- Update clair3 to v0.1.11 for performance improvements
  - --parallel_phase and related options have been removed as intermediate
    phasing no longer needs to be completed in chunks
- longphase is now the default option for intermediate and post-call phasing
- Re-enabled conda profile
- Bumped Aplanat for newer report theme
### Added
- Cumulative read coverage plot added to report

## [v0.2.1]
### Fixed
- Post-call phasing led to mangled VCF.
### Added
- Option to use longphase for post-call phasing.
### Changed
- Disabled conda profile for now.

## [v0.2.0]
### Added
- Singularity profile include in base config.
- Option to use longphase for faster intermediate phasing.
- Improved efficiency of whatshap haplotag.

## [v0.1.2]
### Changed
* Use `bamstats` for alignment statistics instead of `stats_from_bam`.
* BED input now optional (default is all chr1-22,X,Y).
* Update clair3 to v0.1.8 for bugfixes.

## [v0.1.1]
# Fixed
* Correct HTML report filename.

## [v0.1.0]
# Added
* First release.
* Port original Bash script to Nextflow to increase parallelism
* Implement a sharded phasing strategy to improve performance further.
