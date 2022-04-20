## Introduction

This workflow uses [Clair3](https://www.github.com/HKU-BAL/Clair3) for calling small
variants from long reads. Clair3 makes the best of two methods: pileup (for fast
calling of variant candidates in high confidence regions), and full-alignment
(to improve precision of calls of more complex candidates).

The workflow will output a gzipped [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
file containing small variants found in the dataset.


