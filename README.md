# Diploid small variant calling workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for performing diploid variant calling of whole genome data.

## Introduction

This workflow uses [Clair3](https://www.github.com/HKU-BAL/Clair3) for calling small
variants from long reads. Clair3 makes the best of two methods: pileup (for fast
calling of variant candidates in high confidence regions), and full-alignment
(to improve precision of calls of more complex candidates).

The workflow will output a gzipped [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
file containing small variants found in the dataset.


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-human-snp --help
```

to see the workflow options and their descriptions.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software,
it can be downloaded using:

```
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-snp/demo_data.tar.gz
tar -xzvf demo_data.tar.gz
```

The workflow can be run with the demonstration data using:

```
OUTPUT=output
nextflow run epi2me-labs/wf-human-snp \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --bam demo_data/chr6_chr20.bam \
    --bed demo_data/chr6_chr20.bed \
    --ref demo_data/chr6_chr20.fasta \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the 

**Workflow outputs**

The primary outputs of the workflow include:

* a gzipped [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file containing small variants found in the dataset.
* an HTML report document detailing the primary findings of the workflow.
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
