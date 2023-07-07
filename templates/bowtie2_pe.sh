#!/bin/bash

# task.cpus to specify number of threads
# params.genome_fasta_file to specify path reference genome
# name - output file name (provided as arg)

bowtie2 -X2000 --mm -x ${params.bowtie_idx} --threads ${task.cpus} \
    -1 ${fastq1} -2 ${fastq2} \
    | samtools view -Su /dev/stdin \
    | samtools sort - -o ${name}