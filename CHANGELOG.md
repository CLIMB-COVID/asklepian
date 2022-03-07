# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.1.0 2022-03-05
### Changed
* `make_variants_table_v2` is now the default variant table builder
* `gzip --fast` is used to compress the genome table to better balance compression and upload performance

## 3.0.0 2022-02-27
### Changed
* `make_variants_table_v2` uses multiprocessing to parallelize reading the MSA and creating the daily variants table to significantly improve performance

## 2.1.0 2021-08-04
### Added
* `CHANGELOG.md` will now document notable changes to Asklepian
* `environment.yml` provides a conda configuration for the Asklepian environment
### Changed
* `datafunk` has been replaced by `gofasta` to significantly improve performance
* `go.sh` has been updated to better document the global variables required to run, and additionally checks for those variables (exiting if they are not present)
* Genome and variant table generation and uploading subroutines have been moved to `go_genome.sh` and `go_variant.sh` respectively
* `go.sh` now calls the two table subroutines and uses `wait` to run them in parallel instead of serial
* Variables have been tweaked to allow for easier configuration of tests

