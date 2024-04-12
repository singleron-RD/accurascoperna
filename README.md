## Introduction

**singleron-RD/accurascoperna** is a bioinformatics pipeline for processing data from AccuraSCOPE RNA kit: low-throughput single cell 5' + 3' RNA .

In brief, the workflow does the following:

- Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Demultiplexing 3' and 5' reads; convert 5' barcode to 3' barcode; generate STARSolo command-line arguments ([`convert.py`](./bin/convert.py))
- Filter gtf attributes ([`filter_gtf.py`](./bin/filter_gtf.py))
- Generate star genome Index ([`STAR`](https://github.com/alexdobin/STAR/))
- Mapping and quantification of 3' reads. ([`STARSolo`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md))
- Mapping and quantification of 3' read + 5' reads. ([`STARSolo`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md))
- HTML report ([`MultiQC`](http://multiqc.info/))

## Documents

- [Usage](./docs/usage.md)
- [Parameters](./docs/parameters.md)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
