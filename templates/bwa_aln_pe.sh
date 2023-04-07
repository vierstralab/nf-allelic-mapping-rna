#!/bin/bash

# task.cpus to specify number of threads
# params.genome_fasta_file to specify path reference genome
# name - output file name (provided as arg)
# moduleDir - path to directory where pipeline scripts are located
bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
    ${fastq1} \
> pe.reads.rmdup.sorted.remap.fq1.sai

bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
    ${fastq2} \
> pe.reads.rmdup.sorted.remap.fq2.sai

bwa sampe -n 10 -a 750 \
    ${params.genome_fasta_file} \
    pe.reads.rmdup.sorted.remap.fq1.sai pe.reads.rmdup.sorted.remap.fq2.sai \
    pe.reads.rmdup.sorted.remap.fq1.gz pe.reads.rmdup.sorted.remap.fq2.gz \
    | samtools view -b --reference ${params.genome_fasta_file} - \
    > pe.reads.remapped.bam

python3 $moduleDir/bin/filter_reads.py \
    pe.reads.remapped.bam \
    pe.reads.remapped.marked.bam \
    ${params.nuclear_chroms}

samtools sort \
    -@${task.cpus} -l0 pe.reads.remapped.marked.bam \
    | samtools view -b -F 512 - \
    > ${name}