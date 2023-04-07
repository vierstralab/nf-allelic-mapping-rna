#!/bin/bash
# task.cpus to specify number of threads
# params.genome_fasta_file to specify path reference genome
# name - output file name (provided as arg)
# moduleDir - path to directory where pipeline scripts are located

bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
    ${fastq1} > se.reads.rmdup.sorted.remap.fq.sai

bwa samse -n 10 \
    ${params.genome_fasta_file} \
    se.reads.rmdup.sorted.remap.fq.sai \
    se.reads.rmdup.sorted.remap.fq.gz  \
    | samtools view -b --reference ${params.genome_fasta_file} - \
    > se.reads.remapped.bam

python3 $moduleDir/bin/filter_reads.py \
    se.reads.remapped.bam \
    se.reads.remapped.marked.bam \
    ${params.nuclear_chroms}
samtools sort \
        -@${task.cpus} -l0 se.reads.remapped.marked.bam \
    | samtools view -b -F 512 - \
    > ${name}