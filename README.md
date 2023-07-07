# nf-allelic-mapping-rna

This repository hosts a pipeline for allele-specific expression analysis in RNA sequencing data (paired-end). 

### Pipeline Overview
The pipeline consists of the following steps:

- read alignment with [STAR v2.7.9.a](https://github.com/alexdobin/STAR) and simultanious [WASP](https://github.com/bmvdgeijn/WASP) filtering.
- random duplicate removal with [WASP](https://github.com/bmvdgeijn/WASP).
- alignment and WASP filtering/duplicate removal QC.
- gene-level allele-specific expression estimation with [phASER](https://github.com/secastel/phaser).

### Prerequisites
- [nextflow](https://www.nextflow.io/) (22.04.3, dsl2)
- [docker](https://www.docker.com/) / [apptainer](https://apptainer.org) (for containerization)
- [module](https://github.com/cea-hpc/modules/tree/main) with he following modules compiled:
  - STAR/2.7.9a
  - samtools/1.7
  - bcftools/1.12

### Inputs
Main input (--samples_file) is provided as a tsv table with following fields:
 - indiv_id - unique identifier for donor/patient.
 - cell_type - cell type. All libraries with the same combination of cell_type & indiv_id are merged during the alignment step.
 - r1-fastq - ;-separated list of R1 fastqs.
 - r2-fastq - ;-separated list of R2 fastqs (the order should be the same as R1).

Additional inputs include (see params.config):
 - genotype_file - a multi-sample vcf file with genotype information. Should contain genotypes for all patients/donors and names should be identical to indiv_id.
 - star_idx - a directory with STAR-generated genome indexes.
 - genome_fasta_file - genome fasta file.
 - genome_fai_file - gnome fai file.
 - wasp_container - a docker/apptainer container with WASP binaries.
 - wasp_path - path to WASP binaries in the container.
 - phaser_container - a docker/apptainer container with phaser and all dependancies.
 - genes_bed - file with gene coordinates in bed format for phaser.

### Example

```
nextflow run /path/to/nf-allelic-mapping-rna/main.nf -c /path/to/nf-allelic-mapping/nextflow.config -profile Altius --samples_file samples.tsv  --outdir results -resume
```
- nextflow.config is an example of a nextflow config file with slurm/apptainer backend.

### Containers

WASP container:
```
apptainer build wasp.sif Apptainer 
```

phASER container:

from https://github.com/broadinstitute/gtex-pipeline/blob/master/phASER/Dockerfile
```
cd gtex-pipeline/phASER/
docker build --output phaser_docker -t phaser .
```
