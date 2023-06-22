#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
This pipeline counts allele-specific expression for a group of samples.

Input is provided as a table with following fields:
 - indiv_id - unique identifier for donor/patient.
 - cell_type - cell type. All libraries with the same combination of cell_type & indiv_id are merged during the alignment step.
 - r1-fastq - ;-separated list of R1 fastqs.
 - r2-fastq - ;-separated list of R2 fastqs.

Additional inputs include:
 - genotype_file - a multi-sample vcf file with genotype information. Should contain genotypes for all patients/donors and names should be identical to indiv_id.
 - star_idx - a directory with STAR-generated genome indexes.
 - genome_fasta_file - genome fasta file.
 - genome_fai_file - gnome fai file.
 - wasp_container - a docker/apptainer container with WASP binaries.
 - wasp_path - path to WASP binaries in the container.
*/

process generate_vcf {
        // This step is needed because STAR accepts variant information in the form of an uncompressed single-sample vcf.
        // Reheader part and fai-file are needed for GATK compatibility during the count_reads step.
        tag "${indiv_id}:${cell_type}"
        module "bcftools/1.12"
        cpus 10

        input:
        tuple val(indiv_id), val(cell_type), val(fastqs1), val(fastqs2)

        output:
        path("${out_vcf}")

        script:
        out_vcf = "${indiv_id}.${cell_type}.vcf"
        """
        bcftools view -a -s ${indiv_id} ${params.genotype_file} | bcftools view -g het -O v --threads 10 -o intermediate.vcf

        bcftools reheader --fai ${params.genome_fai_file} --threads 10 -o ${out_vcf} intermediate.vcf

        rm intermediate.vcf
        """

}

process align_and_filter {
        /* Notes on STAR parameters:
         - most of the parameters are taken from 'recommended' section in STAR manual.
         - --waspOutputMode SAMtag, --varVCFfile and vW in --outSAMattributes are needed to turn on WASP realignment in STAR.
         - --outFilterScoreMinOverLread 0.3 and --outFilterMatchNminOverLread 0.3 are dealing with excessive filtering of reads and are Alex Dobin-approved.
         - --outFilterMultimapNmax 1 is really important and removes all multi-mappers.
         - --outSAMattrRGline ID:1 SM:sample PL:Illumina PU:x LB:x and --outSAMmapqUnique 60 is for GATK compatibility.
         - BAM is Unsorted in --outSAMtype because it was segfaulting constantly so I moved sorting to a separate step.
        */
        tag "${indiv_id}:${cell_type}"
        module "STAR/2.7.9a:samtools/1.7"
        publishDir params.outdir + '/bams'
        cpus 40
        memory "100 GB"

        input:
        tuple val(indiv_id), val(cell_type), val(fastqs1), val(fastqs2)
        path vcf

        output:
        tuple val(indiv_id), val(cell_type), path("${out_bam}"), path("${out_bam}.bai"), emit: bam
        tuple path('Log.final.out'), path('Aligned.out.bam.flagstat'), path('Aligned.out.filt.bam.flagstat'), emit: qc

        script:
        out_bam = "${indiv_id}.${cell_type}.bam"
        """
        STAR \
                --genomeDir ${params.star_idx} \
                --readFilesIn ${fastqs1} ${fastqs2} \
                --outSAMunmapped Within \
                --outSAMattributes NH HI AS NM MD vW \
                --outFilterMultimapNmax 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --outFilterScoreMinOverLread 0.3 \
                --outFilterMatchNminOverLread 0.3 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --readFilesCommand zcat \
                --runThreadN 40 \
                --genomeLoad LoadAndRemove \
                --limitBAMsortRAM 10000000000 \
                --outSAMtype BAM Unsorted \
                --outSAMheaderHD '@HD' 'VN:1.4' 'SO:coordinate' \
                --outSAMattrRGline ID:1 SM:sample PL:Illumina PU:x LB:x \
                --outSAMmapqUnique 60 \
                --waspOutputMode SAMtag \
                --varVCFfile ${vcf}

        samtools view -h Aligned.out.bam | grep -v "vW:i:[2-7]" | samtools view -b > Aligned.out.filt.bam

        samtools sort -@ 10 -o ${out_bam} Aligned.out.filt.bam

        samtools index ${out_bam}

        samtools flagstat Aligned.out.bam > Aligned.out.bam.flagstat
        samtools flagstat Aligned.out.filt.bam > Aligned.out.filt.bam.flagstat
        """
}        

process rmdup_wasp_style {
        tag "${indiv_id}:${cell_type}"
        container "${params.wasp_container}"
        publishDir params.outdir + '/bams'
        cpus 10

        input:
        tuple val(indiv_id), val(cell_type), path(bam), path(bai)

        output:
        tuple val(indiv_id), val(cell_type), path(out_bam), path("${out_bam}.bai"), emit: bam
        tuple val(indiv_id), val(cell_type), path("${out_bam}.flagstat"), emit: qc

        script:
        out_bam = "${indiv_id}.${cell_type}.rmdup.bam"
        """
        python3 ${params.wasp_path}/mapping/rmdup_pe.py \
                ${bam} rmdup.bam

        samtools sort -@ 10 -o ${out_bam} rmdup.bam
        
        samtools index ${out_bam}

        samtools flagstat ${out_bam} > ${out_bam}.flagstat
        """
}

process count_reads {
        tag "${indiv_id}:${cell_type}"
        module "htslib/1.12:GATK/4.0.1.0"
        publishDir params.outdir + "/read_counts"

        input:
                tuple val(indiv_id), val(cell_type), path(bam), path(bai)
                path vcf

        output:
                tuple val(indiv_id), path(counts)

        script:
        counts = "${indiv_id}.${cell_type}.ase_counts.table"
        """
        bgzip ${vcf}

        tabix ${vcf}.gz

        gatk ASEReadCounter \
                -R ${params.genome_fasta_file} \
                -I ${bam} \
                -V ${vcf}.gz \
                -O ${counts}
        """
}

process compile_qc {
	tag "${indiv_id}:${cell_type}"
	publishDir params.outdir + "/qc"

}

workflow {
        fastqs_grouped_by_lib = Channel
            .fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.indiv_id, row.cell_type, row."r1-fastq", row."r2-fastq") }
            .groupTuple(by: [0,1])
            .map{ it -> tuple(it[0], it[1], it[2].flatten().join(",").replace(";",","), it[3].flatten().join(",").replace(";",",")) }

        vcf = generate_vcf(fastqs_grouped_by_lib)
        bam = align_and_filter(fastqs_grouped_by_lib, vcf)
        rmdup_bam = rmdup_wasp_style(bam.out.bam)
        count_reads(rmdup_bam.out.bam, vcf)
}
