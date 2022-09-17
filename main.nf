#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// TODO move extend metadata process to shared?
// TODO USE DOCKER
params.conda = "/home/sabramov/miniconda3/envs/babachi"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

process generate_h5_tables {
	scratch true
	conda "${params.conda}"
	output:
		path '*.h5'

	script:
	"""
	chroms=("\$(tabix -l ${params.genotype_file})")
	for chrom in \${chroms[@]}; do
		bcftools view -r \${chrom} -Oz ${params.genotype_file} > \${chrom}.vcf.gz
		bcftools index \${chrom}.vcf.gz
	done

	gzip -c ${params.chrom_sizes} > chrom_sizes.txt.gz

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
	conda "${params.conda}"
	publishDir params.outdir + "/remapped"

	cpus 2

	input:
		tuple val(indiv_id), val(ag_number), val(bam_file), val(filtered_sites_file)
		path h5_tables
	
	output:
		tuple val(indiv_id), val(ag_number), val(filtered_sites_file), path("${ag_number}.initial_reads.bed.gz"), file("${ag_number}.initial_reads.bed.gz.tbi"), file("${ag_number}.passing.bam"), file ("${ag_number}.passing.bam.bai")

	script:
	"""
	## split SE from PE reads
	##
	samtools view -O bam -h -F 1 --reference ${params.genome_fasta} ${bam_file} > se.bam
	samtools index se.bam
	n_se=\$(samtools view -c se.bam)

	samtools view -O bam -h -f 1 --reference ${params.genome_fasta} ${bam_file} > pe.bam
	samtools index pe.bam
	n_pe=\$(samtools view -c pe.bam)

	merge_files=""
	remapped_merge_files=""

	## single-ended
	if [[ \${n_se} -gt 0 ]]; then
		
		# an ugly hack to deal with repeated read names on legacy SOLEXA GA1 data
		${moduleDir}/bin/hash_se_reads.py se.bam se.hashed.bam

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

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta} \
			se.reads.rmdup.sorted.remap.fq.gz \
		> se.reads.rmdup.sorted.remap.fq.sai

		bwa samse -n 10 \
			${params.genome_fasta} \
			se.reads.rmdup.sorted.remap.fq.sai \
			se.reads.rmdup.sorted.remap.fq.gz  \
		| samtools view -b --reference ${params.genome_fasta} - \
		> se.reads.remapped.bam

		## step 4 -- mark QC flag
		## Creates filtered bam file se.reads.remapped.marked.bam 

		python3 ${moduleDir}/bin/filter_reads.py \
			se.reads.remapped.bam \
			se.reads.remapped.marked.bam \
			${params.nuclear_chroms}

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
		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta} \
			pe.reads.rmdup.sorted.remap.fq1.gz \
		> pe.reads.rmdup.sorted.remap.fq1.sai

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta} \
			pe.reads.rmdup.sorted.remap.fq2.gz \
		> pe.reads.rmdup.sorted.remap.fq2.sai

		bwa sampe -n 10 -a 750 \
			${params.genome_fasta} \
			pe.reads.rmdup.sorted.remap.fq1.sai pe.reads.rmdup.sorted.remap.fq2.sai \
			pe.reads.rmdup.sorted.remap.fq1.gz pe.reads.rmdup.sorted.remap.fq2.gz \
		| samtools view -b --reference ${params.genome_fasta} - \
		> pe.reads.remapped.bam

		## step 4 -- mark QC flag
		python3 $moduleDir/bin/filter_reads.py \
			pe.reads.remapped.bam \
			pe.reads.remapped.marked.bam \
			${params.nuclear_chroms}

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

	python3 ${moduleDir}/bin/pileup_file.py \
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
	conda "${params.conda}"
	publishDir params.outdir + "/count_reads", mode: 'symlink'

	input:
		tuple val(indiv_id), val(ag_number), val(filtered_sites_file), path(bed_all_reads_file), path(bed_all_reads_file_index), path(bam_passing_file), path(bam_passing_file_index)

	output:
		tuple val(indiv_id), path(name), path("${name}.tbi")

	script:
	name = "${ag_number}.bed.gz"
	"""
	$moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} ${bed_all_reads_file} ${bam_passing_file} | sort-bed - | bgzip -c > ${name}
	
	tabix -p bed ${name}
	"""
}


process merge_by_indiv {
	publishDir params.outdir + "/indiv_merged_files"
	conda "${params.conda}"
	scratch true

	input:
		tuple val(indiv_id), path(bed_files), path(bed_file_index)

	output:
		tuple val(indiv_id), path(name)

	script:
	name = "${indiv_id}.snps.bed"
	"""
	for file in ${bed_files}; do
		python3 $moduleDir/bin/tags_to_babachi_format.py \$file >> ${indiv_id}.snps
	done
	sort -k 1,1 -k2,2n ${indiv_id}.snps > ${name}
	"""
}

workflow waspRealigning {
	take:
		samples_aggregations
	main:
		h5_tables = generate_h5_tables().collect()
		count_reads_files = remap_bamfiles(samples_aggregations, h5_tables) | count_reads
		indiv_merged_count_files = set_key_for_group_tuple(count_reads_files).groupTuple()
		merge_by_indiv(indiv_merged_count_files)
	emit:
		merge_by_indiv.out
}

workflow {
	samples_aggregations = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.indiv_id, row.ag_id, row.bam_file,
		row.filtered_sites_file) }

	waspRealigning(samples_aggregations)
}