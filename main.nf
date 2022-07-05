#!/usr/bin/env nextflow
nextflow.enable.dsl = 1

nuclear_chroms = "${params.genome}.nuclear.txt"
genome_chrom_sizes_file="${params.genome}.chrom_sizes"

genotype_file = "${params.genotyping_output}/genotypes/all.filtered.snps.annotated.vcf.gz"


Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.indiv_id, row.ag_number, row.bamfile, "${params.genotyping_output}/bed_files/${row.filtered_sites_file}", "${params.genotyping_output}/bed_files/${row.filtered_sites_file}.tbi"  ) }
	.tap{ SAMPLES_AGGREGATIONS }

process generate_h5_tables {

	scratch true

	input:
		file vcf_file from file(genotype_file)
		file '*' from file("${genotype_file}.csi")
		file chrom_sizes from file(genome_chrom_sizes_file)

	output:
		file '*.h5' into GENOTYPES_HDF

	script:
	"""
	chroms=("\$(tabix -l ${vcf_file})")
	for chrom in \${chroms[@]}; do
		bcftools view -r \${chrom} -Oz ${vcf_file} > \${chrom}.vcf.gz
		bcftools index \${chrom}.vcf.gz
	done

	gzip -c ${chrom_sizes} > chrom_sizes.txt.gz

	${params.wasp_path}/snp2h5/snp2h5 --chrom chrom_sizes.txt.gz \
		--format vcf \
		--haplotype haplotypes.h5 \
		--snp_index snp_index.h5 \
		--snp_tab snp_tab.h5 \
		chr*.vcf.gz
	"""
}

process remap_bamfiles {
	tag "${indiv_id}:${ag_number}"

	scratch true

	publishDir params.outdir + "/remapped", mode: 'symlink'

	cpus 2

	input:
	set val(indiv_id), val(ag_number), val(bam_file), val(filtered_sites_file), val(filtered_sites_file_index) from SAMPLES_AGGREGATIONS

	file genome from file(params.genome) // doesn't actually make a file
	file '*' from file("${params.genome}.amb")
  	file '*' from file("${params.genome}.ann")
  	file '*' from file("${params.genome}.bwt")
  	file '*' from file("${params.genome}.fai")
  	file '*' from file("${params.genome}.pac")
  	file '*' from file("${params.genome}.sa")
	
	file nuclear_chroms from file(nuclear_chroms)
	file '*' from GENOTYPES_HDF.collect()
	
	output:
	set val(indiv_id), val(ag_number), val(filtered_sites_file), val(filtered_sites_file_index), path("${ag_number}.initial_reads.bed.gz"), file("${ag_number}.initial_reads.bed.gz.tbi"), file("${ag_number}.passing.bam"), file ("${ag_number}.passing.bam.bai") into REMAPPED_READS

	script:
	"""
	## split SE from PE reads
	##
	samtools view -O bam -h -F 1 ${bam_file} > se.bam
	samtools index se.bam
	n_se=\$(samtools view -c se.bam)

	samtools view -O bam -h -f 1 ${bam_file} > pe.bam
	samtools index pe.bam
	n_pe=\$(samtools view -c pe.bam)

	merge_files=""
	remapped_merge_files=""

	## single-ended
	if [[ \${n_se} -gt 0 ]]; then
		
		# an ugly hack to deal with repeated read names on legacy SOLEXA GA1 data
		$baseDir/bin/hash_se_reads.py se.bam se.hashed.bam

		## step 1 -- remove duplicates
		##
		python3 ${params.wasp_path}/mapping/rmdup.py \
			se.hashed.bam  se.reads.rmdup.bam

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} \
			-o se.reads.rmdup.sorted.bam \
			-O bam \
			se.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		### Creates 3 following files:
		### se.reads.rmdup.sorted.to.remap.bam (reads to remap)
		### se.reads.rmdup.sorted.keep.bam (reads to keep)
		### se.reads.rmdup.sorted.remap.fq.gz (fastq file containing the reads with flipped alleles to remap)
		python3 ${params.wasp_path}/mapping/find_intersecting_snps.py \
			--is_sorted \
			--output_dir \${PWD} \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			se.reads.rmdup.sorted.bam

		## step 3 -- re-align reads
		##

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${genome} \
			se.reads.rmdup.sorted.remap.fq.gz \
		> se.reads.rmdup.sorted.remap.fq.sai

		bwa samse -n 10 \
			${genome} \
			se.reads.rmdup.sorted.remap.fq.sai \
			se.reads.rmdup.sorted.remap.fq.gz  \
		| samtools view -b -t ${genome}.fai - \
		> se.reads.remapped.bam

		## step 4 -- mark QC flag
		## Creates filtered bam file se.reads.remapped.marked.bam 

		python3 $baseDir/bin/filter_reads.py \
			se.reads.remapped.bam \
			se.reads.remapped.marked.bam \
			${nuclear_chroms}

		## step 5 -- filter reads

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} -l0 se.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> se.reads.remapped.marked.filtered.bam

		python3 ${params.wasp_path}/mapping/filter_remapped_reads.py \
			se.reads.rmdup.sorted.to.remap.bam \
			se.reads.remapped.marked.filtered.bam \
			se.reads.remapped.result.bam

		merge_files="\${merge_files} se.reads.rmdup.sorted.bam"
		remapped_merge_files="\${remapped_merge_files} se.reads.remapped.result.bam se.reads.rmdup.sorted.keep.bam"
	fi

	## paired-end
	if [[ \${n_pe} -gt 0 ]]; then
		
		## step 1 -- remove duplicates
		##
		python3 ${params.wasp_path}/mapping/rmdup_pe.py \
			pe.bam pe.reads.rmdup.bam

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} \
			-o pe.reads.rmdup.sorted.bam \
			-O bam \
			pe.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		##
		python3 ${params.wasp_path}/mapping/find_intersecting_snps.py \
			--is_paired_end \
			--is_sorted \
			--output_dir \${PWD} \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			pe.reads.rmdup.sorted.bam

		## step 3 -- re-align reads
		##

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${genome} \
			pe.reads.rmdup.sorted.remap.fq1.gz \
		> pe.reads.rmdup.sorted.remap.fq1.sai

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${genome} \
			pe.reads.rmdup.sorted.remap.fq2.gz \
		> pe.reads.rmdup.sorted.remap.fq2.sai

		bwa sampe -n 10 -a 750 \
			${genome} \
			pe.reads.rmdup.sorted.remap.fq1.sai pe.reads.rmdup.sorted.remap.fq2.sai \
			pe.reads.rmdup.sorted.remap.fq1.gz pe.reads.rmdup.sorted.remap.fq2.gz \
		| samtools view -b -t ${genome}.fai - \
		> pe.reads.remapped.bam

		## step 4 -- mark QC flag
		##

		python3 $baseDir/bin/filter_reads.py \
			pe.reads.remapped.bam \
			pe.reads.remapped.marked.bam \
			${nuclear_chroms}

		## step 5 -- filter reads

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} -l0 pe.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> pe.reads.remapped.marked.filtered.bam

		python3 ${params.wasp_path}/mapping/filter_remapped_reads.py \
			pe.reads.rmdup.sorted.to.remap.bam \
			pe.reads.remapped.marked.filtered.bam \
			pe.reads.remapped.result.bam

		merge_files="\${merge_files} pe.reads.rmdup.sorted.bam"
		remapped_merge_files="\${remapped_merge_files} pe.reads.remapped.result.bam"
	fi



	## step 6 -- merge back reads
	samtools merge -f \
		reads.rmdup.bam \
		\${merge_files}

	samtools sort \
		-m ${task.memory.toMega()}M \
		-@${task.cpus} \
		-o reads.rmdup.sorted.bam  \
		reads.rmdup.bam
	
	samtools index reads.rmdup.sorted.bam

	python3 $baseDir/bin/pileup_file.py \
		 ${filtered_sites_file} reads.rmdup.sorted.bam \
		  | sort-bed - | bgzip -c > ${ag_number}.initial_reads.bed.gz
	
	tabix -p bed ${ag_number}.initial_reads.bed.gz
	# todo: merge dedupped se and pe reads

	samtools merge -f reads.passing.bam \
		\${remapped_merge_files}

	samtools sort \
		-m ${task.memory.toMega()}M \
		-@${task.cpus} \
		-o reads.passing.sorted.bam  \
		reads.passing.bam 

	mv reads.passing.sorted.bam ${ag_number}.passing.bam
	samtools index ${ag_number}.passing.bam
	"""
}

process count_reads {
	tag "${indiv_id}:${ag_number}"

	publishDir params.outdir + "/count_reads", mode: 'symlink'

	input:
	set val(indiv_id), val(ag_number), val(filtered_sites_file), val(filtered_sites_file_index), file(bed_all_reads_file), file(bed_all_reads_file_index), file(bam_passing_file), file(bam_passing_file_index) from REMAPPED_READS

	output:
	set val(indiv_id), file(name), file("${name}.tbi") into COUNT_READS_FILES

	script:
	name = "${ag_number}.bed.gz"
	"""
	$baseDir/bin/count_tags_pileup.py \
		${filtered_sites_file} ${bed_all_reads_file} ${bam_passing_file} \
		 | sort-bed - | bgzip -c > ${name}
	
	tabix -p bed ${name}
	"""
}

INDIV_MERGED_COUNT_FILES = COUNT_READS_FILES.groupTuple().map{ it -> tuple(it[0], it[1].join(" "), it[2].join(" ")) }

process merge_by_indiv {
	publishDir params.outdir + "/indiv_merged_files", mode: 'symlink'

	input:
	tuple val(indiv_id), path(bed_files) from INDIV_MERGED_COUNT_FILES

	output:
	tuple val(indiv_id), file(name)

	script:
	name = "${indiv_id}.snps.bed"
	"""
	for file in ${bed_files}; do
		python3 $baseDir/bin/tags_to_babachi_format.py \$file >> ${indiv_id}.snps
	done
	sort -k 1,1 -k2,2n ${indiv_id}.snps > ${name}
	rm ${indiv_id}.snps
	"""
}
