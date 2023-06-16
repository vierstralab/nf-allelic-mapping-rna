#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.star_idx = "/net/seq/data/genomes/human/GRCh38/noalts-sequins/STARgenome-gencode-v25/"
params.genotype_file = "/home/ezotova/230227_regulotyping/data/genotyping_panel/chroms1-22.phaseI+II+stim.annotated.ancestral.dedup.gene_regions_only.vcf.gz"

process align_and_filter {
	tag "${indiv_id}:${cell_type}"
	module "bedtools/2.25.0:STAR/2.7.10a:samtools/1.3"
	publishDir params.outdir + '/bams'
	cpus 40

	input:
	tuple val(indiv_id), val(cell_type), val(fastqs1), val(fastqs2)

	output:
	tuple val(indiv_id), val(cell_type), path('*.bam'), path('*.bai'), emit: tuple_with_bam
	path("${indiv_id}.vcf"), emit: vcf 

	script:
	"""
	bcftools view -a -s ${indiv_id} ${params.genotype_file} | bcftools view -g het -O v --threads 10 -o ${indiv_id}.vcf

	STAR \
		--genomeDir ${params.star_idx} \
		--readFilesIn ${fastqs1} ${fastqs2} \
		--outSAMunmapped Within --outFilterType BySJout \
		--outSAMattributes NH HI AS NM MD vW \
		--outFilterMultimapNmax 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--sjdbScore 1 \
		--readFilesCommand zcat \
		--runThreadN 40 \
		--genomeLoad LoadAndRemove \
		--limitBAMsortRAM 10000000000 \
		--outSAMtype BAM Unsorted \
		--outSAMheaderHD '@HD' 'VN:1.4' 'SO:coordinate' \
		--waspOutputMode SAMtag \
		--varVCFfile ${indiv_id}.vcf

	samtools view -h Aligned.out.bam | grep -v "vW:i:[2-7]" | samtools view -b  > Aligned.out.filt.bam

	samtools sort  -@ 40  -o ${indiv_id}.bam  Aligned.out.filt.bam

	samtools index ${indiv_id}.bam
	"""
}	

process rmdup_wasp_style {
	tag "${indiv_id}:${cell_type}"
	container "${params.wasp_container}"
	publishDir params.outdir + '/bams'
	cpus 2

	input:
	tuple val(indiv_id), val(cell_type), path(bam), path(bai)

	output:
	tuple val(indiv_id), val(cell_type), path(out_bam), path("${out_bam}.bai")

	script:
	name = "${indiv_id}.${cell_type}.rmdup"
	out_bam = "${name}.bam"
	"""
	python3 ${params.wasp_path}/mapping/rmdup_pe.py \
		${bam} ${out_bam}
	
	samtools index ${out_bam}
	"""
}

process count_reads {
        tag "${indiv_id}:${cell_type}"
	module "GATK/4.0.1.0"
        publishDir params.outdir + "/read_counts"

        input:
                tuple val(indiv_id), val(cell_type), path(bam), path(bai)
		tuple path(vcf)

        output:
                tuple val(indiv_id), path(counts)

        script:
        counts = "${indiv_id}.${cell_type}.ase_counts.table"
        """
        gatk ASEReadCounter \
                -R ${params.genome_fasta_file} \
                -I ${bam} \
                -V ${vcf} \
                -O ${counts}
        """
}

workflow {
	fastqs_grouped_by_lib = Channel
	    .fromPath(params.samples_file)
	    .splitCsv(header:true, sep:'\t')
	    .map{ row -> tuple(row.indiv_id, row.cell_type, row."r1-fastq", row."r2-fastq") }
	    .groupTuple(by: [0,1])
	    .map{ it -> tuple(it[0], it[1], it[2].flatten().join(",").replace(";",","), it[3].flatten().join(",").replace(";",",")) }

	vcf_and_bam = align_and_filter(fastqs_grouped_by_lib)
	rmdup_bam = rmdup_wasp_style(vcf_and_bam.tuple_with_bam)
	count_reads(rmdup_bam, vcf_and_bam.vcf)
}
