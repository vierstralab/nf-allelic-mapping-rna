#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

def get_container(file_name) {
  parent = file(file_name).parent
  container = "-v ${parent}:${parent}"
}

conda = '/home/sabramov/miniconda3/envs/babachi/envs/allelic-mapping'
//wasp_path = '/opt/WASP'
wasp_path = '/home/sabramov/projects/ENCODE4/WASP'


process filter_variants {
	tag "${indiv_id}"
	conda conda
	publishDir "${params.outdir}/bed_files"

	input:
		tuple val(ag_id), val(indiv_id)

	output:
		tuple val(ag_id), path(outname), path("${outname}.tbi")

	script:
	outname = "${indiv_id}:${ag_id}.bed.gz"
	"""
	bcftools query \
		-s ${indiv_id} \
		-i'GT="alt"' \
		-f'%CHROM\\t%POS0\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/MAF\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
		${params.genotype_file} \
	| awk -v OFS="\\t" \
		-v min_GQ=${params.min_GQ} -v min_AD=${params.min_AD} -v min_DP=${params.min_DP}\
		'\$9<min_GQ { next; } \$10<min_DP { next; }\
			(\$8=="0/1" || \$8=="1/0" || \$8=="0|1" || \$8=="1|0") \
			&& (\$11<min_AD || \$12<min_AD) { next; } \
			{ print; }' \
	| sort-bed - \
	| { grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn || true; } \
	| bgzip -c > ${outname}
	tabix -f -p bed "${outname}"
	"""
}


process generate_h5_tables {
	scratch true
	// container "${params.container}"
	// containerOptions "${get_container(params.genotype_file)} ${get_container(params.chrom_sizes)}"
	conda conda

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

	${wasp_path}/snp2h5/snp2h5 --chrom chrom_sizes.txt.gz \
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
	conda conda
	// container "${params.container}"
	// containerOptions "${get_container(params.genome_fasta)} ${get_container(params.nuclear_chroms)}"
	publishDir params.outdir + "/remapped", pattern: "${ag_number}.passing.bam*"

	cpus 2

	input:
		tuple val(ag_number), val(indiv_id), val(bam_file), path(filtered_sites_file), path(filtered_sites_file_index)
		path h5_tables
	
	output:
		tuple val(indiv_id), val(ag_number), path(filtered_sites_file), path(filtered_sites_file_index), path("${ag_number}.passing.bam"), path("${ag_number}.passing.bam.bai")

	script:
	mem = Math.round(task.memory.toMega() / task.cpus * 0.95)
	"""
	## split SE from PE reads
	##
	samtools view -O bam -h -F 1 --reference ${params.genome_fasta} ${bam_file} > se.bam
	samtools index se.bam
	n_se=\$(samtools view -c se.bam)

	samtools view -O bam -h -f 1 --reference ${params.genome_fasta} ${bam_file} > pe.bam
	samtools index pe.bam
	n_pe=\$(samtools view -c pe.bam)

	remapped_merge_files=""

	## single-ended
	if [[ \${n_se} -gt 0 ]]; then
		
		# an ugly hack to deal with repeated read names on legacy SOLEXA GA1 data
		$moduleDir/bin/hash_se_reads.py se.bam se.hashed.bam

		## step 1 -- remove duplicates
		##
		python3 ${wasp_path}/mapping/rmdup.py \
			se.hashed.bam  se.reads.rmdup.bam

		samtools sort \
			-m ${mem}M \
			-@${task.cpus} \
			-o se.reads.rmdup.sorted.bam \
			-O bam \
			se.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		### Creates 3 following files:
		### se.reads.rmdup.sorted.to.remap.bam (reads to remap)
		### se.reads.rmdup.sorted.keep.bam (reads to keep)
		### se.reads.rmdup.sorted.remap.fq.gz (fastq file containing the reads with flipped alleles to remap)
		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
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

		python3 $moduleDir/bin/filter_reads.py \
			se.reads.remapped.bam \
			se.reads.remapped.marked.bam \
			${params.nuclear_chroms}

		## step 5 -- filter reads

		samtools sort \
			-m ${mem}M \
			-@${task.cpus} -l0 se.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> se.reads.remapped.marked.filtered.bam

		python3 ${wasp_path}/mapping/filter_remapped_reads.py \
			se.reads.rmdup.sorted.to.remap.bam \
			se.reads.remapped.marked.filtered.bam \
			se.reads.remapped.result.bam
		
		remapped_merge_files="\${remapped_merge_files} se.reads.remapped.result.bam"
	fi

	## paired-end
	if [[ \${n_pe} -gt 0 ]]; then
		
		## step 1 -- remove duplicates
		##
		python3 ${wasp_path}/mapping/rmdup_pe.py \
			pe.bam pe.reads.rmdup.bam

		samtools sort \
			-m ${mem}M \
			-@${task.cpus} \
			-o pe.reads.rmdup.sorted.bam \
			-O bam \
			pe.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		##
		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
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
			-m ${mem}M \
			-@${task.cpus} -l0 pe.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> pe.reads.remapped.marked.filtered.bam

		python3 ${wasp_path}/mapping/filter_remapped_reads.py \
			pe.reads.rmdup.sorted.to.remap.bam \
			pe.reads.remapped.marked.filtered.bam \
			pe.reads.remapped.result.bam
		remapped_merge_files="\${remapped_merge_files} pe.reads.remapped.result.bam"
	fi



	## step 6 -- merge back reads

	samtools merge -f reads.passing.bam \
		\${remapped_merge_files}

	samtools sort \
		-m ${mem}M \
		-@${task.cpus} \
		-o reads.passing.sorted.bam  \
		reads.passing.bam 

	mv reads.passing.sorted.bam ${ag_number}.passing.bam
	samtools index ${ag_number}.passing.bam
	"""
}

process count_reads {
	tag "${indiv_id}:${ag_number}"
	// container "${params.container}"
	conda conda
	publishDir params.outdir + "/count_reads"

	input:
		tuple val(indiv_id), val(ag_number), path(filtered_sites_file), path(filtered_sites_file_index), path(bam_passing_file), path(bam_passing_file_index)

	output:
		tuple val(indiv_id), path(name), path("${name}.tbi")

	script:
	name = "${ag_number}.bed.gz"
	"""
	$moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} ${bam_passing_file} | sort-bed - | bgzip -c > ${name}
	
	tabix -p bed ${name}
	"""
}


process merge_by_indiv {
	publishDir params.outdir + "/indiv_merged_files.maf_filter"
	tag "${indiv_id}"
	// container "${params.container}"
	conda conda
	scratch true

	input:
		tuple val(indiv_id), path(bed_files), path(bed_file_index)

	output:
		tuple val(indiv_id), path(name)

	script:
	name = "${indiv_id}.snps.bed"
	"""
	for file in ${bed_files}
	do
		python3 $moduleDir/bin/tags_to_babachi_format.py \$file >> ${indiv_id}.snps
	done
	sort -k 1,1 -k2,2n ${indiv_id}.snps > ${name}
	"""
}

workflow waspRealigning {
	take:
		samples_aggregations
	main:
		sagr = samples_aggregations.map(it -> tuple(it[1], it[0], it[2]))
		h5_tables = generate_h5_tables().collect()
		snps_sites = filter_variants(sagr.map(it -> tuple(it[0], it[1])))
		samples = sagr.join(snps_sites, by: 0)
		count_reads_files = remap_bamfiles(samples, h5_tables) | count_reads
		indiv_merged_count_files = count_reads_files.groupTuple()
		merge_by_indiv(indiv_merged_count_files)
	emit:
		merge_by_indiv.out
}

workflow test {
	base_path = '/net/seq/data2/projects/sabramov/ENCODE4/wasp-realigning/output'
	val(indiv_id), path(bed_files), path(bed_file_index)
	samples_aggregations = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple(row.indiv_id,
			file("${base_path}/count_reads/${row.ag_id}.bed.gz"),
			file("${base_path}/count_reads/${row.ag_id}.bed.gz.tbi")))
		.unique { it[1] }
	
	merge_by_indiv(samples_aggregations.groupTuple())
	
}

workflow {
	samples_aggregations = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple(row.indiv_id, row.ag_id, row.bam_file)).unique { it[1] }
	waspRealigning(set_key_for_group_tuple(samples_aggregations))
}