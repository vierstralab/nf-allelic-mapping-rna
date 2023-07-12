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
 - phaser_container - a docker/apptainer container with phaser and all dependancies.
 - genes_bed - file with gene coordinates in bed format for phaser.
*/

process generate_vcf {
        // This step is needed because STAR accepts variant information in the form of an uncompressed single-sample vcf.
        // Reheader part and fai-file are needed for GATK compatibility during the count_reads step.
        tag "${indiv_id}"
        module "bcftools/1.12"
        cpus 10

        input:
        val indiv_id

        output:
        tuple val("${indiv_id}"), path("${out_vcf}")

        script:
        out_vcf = "${indiv_id}.vcf"
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
        tuple val(indiv_id), val(cell_type), val(fastqs1), val(fastqs2), path(vcf)

        output:
        tuple val(indiv_id), val(cell_type), path("${out_bam}"), path("${out_bam}.bai"), emit: bam
        tuple val(indiv_id), val(cell_type), path('Log.final.out'), path('Aligned.out.bam.flagstat'), path('Aligned.out.filt.bam.flagstat'), emit: qc

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
        tuple val(indiv_id), val(cell_type), path("${out_bam}.flagstat"), path('.command.log'), emit: qc

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
        container "${params.phaser_container}"
        cpus 20

        publishDir params.outdir + "/read_counts"

        input:
                tuple val(indiv_id), val(cell_type), path(bam), path(bai), path(vcf)

        output:
                path("${indiv_id}.${cell_type}.allel*")
                path("${indiv_id}.${cell_type}.gene_ae_counts.txt")

        script:
        """
        bgzip ${vcf}

        tabix ${vcf}.gz

        python2.7 /opt/phaser/phaser/phaser.py --vcf ${vcf}.gz --bam ${bam} --paired_end 1 --map 60 --baseq 10 --sample ${indiv_id} --o ${indiv_id}.${cell_type} --threads 20

        python2.7 /opt/phaser/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts ${indiv_id}.${cell_type}.haplotypic_counts.txt --features ${params.genes_bed} --o ${indiv_id}.${cell_type}.gene_ae_counts.txt
        """
}

process compile_qc {
        /*
        This process parces various logs with awk and collects the metrics into a sindle-line tsv table with meaningful header.
        Flagstat files produced by different versions of samtools have different formats.
        */
        tag "${indiv_id}:${cell_type}"
        module "samtools/1.7"
        publishDir params.outdir + "/qc"
        input:
                tuple val(indiv_id), val(cell_type), path(rmdup_flagstat), path(rmdup_log), path(STAR_log), path(raw_flagstat), path(filt_flagstat), path(rmdup_bam), path(rmdup_bai)

        output:
                path("${out}")

        shell:
        out="${indiv_id}.${cell_type}.qc.tsv"
        '''
        echo -e "\
indiv_id\tcell_type\t\
raw_flagstat_total_reads\t\
raw_flagstat_%_mapped_reads\t\
raw_star_%_removed_multimappers\t\
raw_star_%_removed_by_length_score\t\
wasp_filt_flagstat_total_reads\t\
wasp_filt_flagstat_%_mapped_reads\t\
wasp_filt_and_dedup_flagstat_total_reads\t\
wasp_filt_and_dedup_flagstat_%_mapped_reads\t\
wasp_filt_and_dedup_log_removed_unaligned_readpairs\t\
wasp_filt_and_dedup_log_removed_duplicated_readpairs\t\
wasp_filt_and_dedup_log_kept_readpairs\t\
wasp_filt_and_dedup_sequins_mapped_reads
" > !{out}

        echo -e "\
!{indiv_id}\t!{cell_type}\t\
$(awk 'NR==1 {print $1}' !{raw_flagstat})\t\
$(awk 'NR==5 {split($5, a, "%"); sub(/[()]/, "", a[1]); print a[1]}' !{raw_flagstat})\t\
$(grep "% of reads mapped to too many loci" !{STAR_log} | awk '{sub(/%/,"",$10); print $10}')\t\
$(grep "% of reads unmapped: too short" !{STAR_log} | awk '{sub(/%/,"",$8); print $8}')\t\
$(awk 'NR==1 {print $1}' !{filt_flagstat})\t\
$(awk 'NR==5 {split($5, a, "%"); sub(/[()]/, "", a[1]); print a[1]}' !{filt_flagstat})\t\
$(awk 'NR==1 {print $1}' !{rmdup_flagstat})\t\
$(awk 'NR==7 {split($5, a, "%"); sub(/[()]/, "", a[1]); print a[1]}' !{rmdup_flagstat})\t\
$(grep "  unmapped:" !{rmdup_log} | awk '{print $2}')\t\
$(grep "  duplicate pairs:" !{rmdup_log} | awk '{print $3}')\t\
$(grep "  pairs:" !{rmdup_log} | awk '{print $2}')\t\
$(samtools view !{rmdup_bam} | grep -c chrIS)
" >> !{out}
        '''
}


workflow {
        fastqs_grouped_by_lib = Channel
            .fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.indiv_id, row.cell_type, row."r1-fastq", row."r2-fastq") }
            .groupTuple(by: [0,1])
            .map{ it -> tuple(it[0], it[1], it[2].flatten().join(",").replace(";",","), it[3].flatten().join(",").replace(";",",")) }

        vcf = fastqs_grouped_by_lib
		.map{ it -> it[0] }
		.unique()
		| generate_vcf

        bam = vcf
		.cross(fastqs_grouped_by_lib)
		.map{ it -> it[1..0].flatten()[0..3,5] }  // indiv_id, cell_type, r1-fastq, r2-fastq, vcf
		| align_and_filter

        rmdup_bam = rmdup_wasp_style(bam.bam)

        vcf.cross(rmdup_bam.bam)
		.map{ it -> it[1..0].flatten()[0..3,5] } // indiv_id, cell_type, bam, bai, vcf
		| count_reads

	rmdup_bam.qc
		.join( bam.qc, by: [0,1] )
		.join( rmdup_bam.bam, by: [0,1] )
		| compile_qc
}
