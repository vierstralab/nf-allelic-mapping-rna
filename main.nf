#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

def get_container(file_name) {
  parent = file(file_name).parent
  old_parent = file(file_name).toRealPath().parent
  container = "--bind ${parent},${old_parent}"
}

wasp_path = '/opt/WASP'


process make_iupac_genome {
	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)} ${get_container(params.genotype_file)}" 
	publishDir "${params.outdir}/alt_genome"
	
	output:
		tuple path("${name}"), path("${name}.fai")

	script:
	name = "iupac.genome.fa"
    """
    python3 $moduleDir/bin/nonref_genome.py ${params.genome_fasta_file} ${name} ${params.genotype_file}
    """
}


process filter_variants {
	tag "${indiv_id}:${ag_id}"
	container "${params.container}"
	containerOptions "${get_container(params.genotype_file)}" 
	publishDir "${params.outdir}/target_variants"

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
		-f"%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t%INFO/AAF\t%INFO/RAF\t[%GT\t%GQ\t%DP\t%AD{0}\t%AD{1}]\n" \
		${params.genotype_file} \
	| awk -v OFS="\\t" \
		-v min_GQ=${params.min_GQ} -v min_AD=${params.min_AD} -v min_DP=${params.min_DP}\
		'\$10<min_GQ { next; } \$11<min_DP { next; }\
			(\$9=="0/1" || \$9=="1/0" || \$9=="0|1" || \$9=="1|0") \
			&& (\$12<min_AD || \$13<min_AD) { next; } \
			{ print; }' \
	| sort-bed - \
	| { grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn || true; } \
	| bgzip -c > ${outname}
	tabix -f -p bed "${outname}"
	"""
}


process generate_h5_tables {
	scratch true
	container "${params.container}"
	containerOptions "${get_container(params.genotype_file)} ${get_container(params.chrom_sizes)}"

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
	scratch "$workDir"
	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)} ${get_container(params.nuclear_chroms)}"
	cpus 2

	input:
		tuple val(ag_number), val(indiv_id), path(bam_file), path(bam_index_file), path(filtered_sites_file), path(filtered_sites_file_index)
		path h5_tables
	
	output:
		tuple val(indiv_id), val(ag_number), path(filtered_sites_file), path(filtered_sites_file_index), path("${ag_number}.passing.bam"), path("${ag_number}.passing.bam.bai"), path("${ag_number}.coverage.bed.gz"), path("${ag_number}.coverage.bed.gz.tbi")

	script:
	mem = Math.round(task.memory.toMega() / task.cpus * 0.95)
	"""
	## split SE from PE reads
	##
	samtools view -O bam -h -F 1 --reference ${params.genome_fasta_file} ${bam_file} > se.bam
	samtools index se.bam
	n_se=\$(samtools view -c se.bam)

	samtools view -O bam -h -f 1 --reference ${params.genome_fasta_file} ${bam_file} > pe.bam
	samtools index pe.bam
	n_pe=\$(samtools view -c pe.bam)

	remapped_merge_files=""
	rmdup_original_files=""

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
		
		rmdup_original_files="\${rmdup_original_files} se.reads.rmdup.sorted.bam"

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

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
			se.reads.rmdup.sorted.remap.fq.gz \
		> se.reads.rmdup.sorted.remap.fq.sai

		bwa samse -n 10 \
			${params.genome_fasta_file} \
			se.reads.rmdup.sorted.remap.fq.sai \
			se.reads.rmdup.sorted.remap.fq.gz  \
		| samtools view -b --reference ${params.genome_fasta_file} - \
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
		
		rmdup_original_files="\${rmdup_original_files} pe.reads.rmdup.sorted.bam"

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
		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
			pe.reads.rmdup.sorted.remap.fq1.gz \
		> pe.reads.rmdup.sorted.remap.fq1.sai

		bwa aln -Y -l 32 -n 0.04 -t ${task.cpus} ${params.genome_fasta_file} \
			pe.reads.rmdup.sorted.remap.fq2.gz \
		> pe.reads.rmdup.sorted.remap.fq2.sai

		bwa sampe -n 10 -a 750 \
			${params.genome_fasta_file} \
			pe.reads.rmdup.sorted.remap.fq1.sai pe.reads.rmdup.sorted.remap.fq2.sai \
			pe.reads.rmdup.sorted.remap.fq1.gz pe.reads.rmdup.sorted.remap.fq2.gz \
		| samtools view -b --reference ${params.genome_fasta_file} - \
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



	## step 6 -- merge se and pe reads
	if [ "`echo \${remapped_merge_files} | wc -w`" -ge 2 ]; then
		samtools merge -f reads.passing.bam \
			\${remapped_merge_files}
	else
		mv \${remapped_merge_files} reads.passing.bam
	fi
		
	samtools sort \
		-m ${mem}M \
		-@${task.cpus} \
		-o ${ag_number}.passing.bam  \
		reads.passing.bam
	samtools index ${ag_number}.passing.bam


	###########################
	if [ "`echo \${rmdup_original_files} | wc -w`" -ge 2 ]; then
		samtools merge -f rmdup_original_files \
			\${rmdup_original_files}

		samtools sort \
		-m ${mem}M \
		-@${task.cpus} \
		-o reads.original.sorted.rmdup.bam  \
		reads.rmdup.original.bam
	else
		mv \${rmdup_original_files} reads.original.sorted.rmdup.bam
	fi

	samtools index reads.original.sorted.rmdup.bam

	python3 $moduleDir/bin/count_tags_pileup.py ${filtered_sites_file} \
		 reads.original.sorted.rmdup.bam \
		 --only_coverage | bgzip -c > ${ag_number}.coverage.bed.gz
	tabix ${ag_number}.coverage.bed.gz
	"""
}

process count_reads {
	tag "${indiv_id}:${ag_number}"
	container "${params.container}"
	publishDir params.outdir + "/count_reads"

	input:
		tuple val(indiv_id), val(ag_number), path(filtered_sites_file), path(filtered_sites_file_index), path(bam_passing_file), path(bam_passing_file_index), path(rmdup_counts), path(rmdup_counts_index)

	output:
		tuple val(indiv_id), path(name), path("${name}.tbi")

	script:
	name = "${ag_number}.bed.gz"
	"""
	python3 $moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} ${bam_passing_file} --original_dedup_cover ${rmdup_counts} | sort-bed - | bgzip -c > ${name}
	tabix ${name}
	"""
}

process merge_by_indiv {
	publishDir "${params.outdir}/indiv_merged_files"
	tag "${indiv_id}"
	container "${params.container}"
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
	echo "#chr\tstart\tend\tID\tref\talt\tref_counts\talt_counts\tsample_id\tAAF\tRAF\tFMR" > ${name}
	sort -k 1,1 -k2,2n ${indiv_id}.snps >> ${name}
	"""
}

process add_snp_files_to_meta {
	publishDir "${params.outdir}"
	container "${params.container}"
	containerOptions "${get_container(params.samples_file)}"

	output:
		path name

	script:
	name = "meta+sample_ids.tsv"
	"""
	python3 $moduleDir/bin/add_meta.py ${params.samples_file} ${name} ${launchDir}/${params.outdir}/indiv_merged_files
	"""
}


workflow waspRealigning {
	take:
		samples_aggregations
	main:
		sagr = samples_aggregations.map(it -> tuple(it[1], it[0], it[2], it[3]))
		h5_tables = generate_h5_tables().collect()
		snps_sites = filter_variants(sagr.map(it -> tuple(it[0], it[1])))
		samples = sagr.join(snps_sites, by: 0)
		count_reads_files = remap_bamfiles(samples, h5_tables) | count_reads
		out = merge_by_indiv(count_reads_files.groupTuple())
	emit:
		out
}


workflow makeAltGenome {
	make_iupac_genome()
}

workflow {
	samples_aggregations = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple(row.indiv_id, row.ag_id, file(row.bam_file), file("${row.bam_file}.crai")))
		.unique { it[1] }
	indivs_count = samples_aggregations.map(it -> it[0]).unique().count().view {
		it -> """There are ${it} unique INDIV_IDs in the ${params.samples_file}. Please, check that they correspond to IDs in ${params.genotype_file}"""
	}
	waspRealigning(set_key_for_group_tuple(samples_aggregations))
	add_snp_files_to_meta() 
}