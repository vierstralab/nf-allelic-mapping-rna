#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ADD more aligners in the process below if needed. See example in the script body
process align_reads {
	container "${params.container}"
	// All external files (passed as params) should be wrapped 
	// in "get_container" function and added to containerOptions
	containerOptions "${get_container(params.genome_fasta_file)} ${get_container(params.nuclear_chroms)} ${params.aligner == 'bowtie-chip' ? get_container(params.bowtie_idx) : ''}"
	scratch true
	tag "${ag_number}:${r_tag}"
	cpus 3
	// fastq1 is the same as fastq2 if r_tag == "se"
	input:
		tuple val(ag_number), val(r_tag), path(fastq1), path(fastq2)
		// "ag_number" - ID of initial bam file that was splitted
		// into pe and se reads.

		// "r_tag" - indicates which reads were extracted from the bam file.
		// Can be either "se" or "pe"

		// fastq1, fastq2 reads from the bam file. 
		// Use any of them for se reads aligning

	output:
		tuple val(ag_number), val(r_tag), path(name)

	script:
	// Name of the output bam file
	name = "${ag_number}.${r_tag}.realigned.bam"
	switch(params.aligner) {
		case 'bwa-altius-dnase':
			// PE reads alignment
			if (r_tag == 'pe') {
				"""
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
				""" 
			} else {
				// SE reads alignment
				"""
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
				"""
			};
			break;
		case "bowtie-chip":
			if (r_tag == 'pe') {
				"""
				bowtie2 -X2000 --mm -x ${params.bowtie_idx} --threads ${task.cpus} \
   		 			-1 ${fastq1} -2 ${fastq2} \
		 			| samtools view -Su /dev/stdin \
					| samtools sort - -o ${name}
				"""
			} else {
				"""
				bowtie2 --mm -x ${params.bowtie_idx} --threads ${task.cpus} \
					-U <(zcat -f ${fastq1}) \
					| samtools view -Su /dev/stdin \
					| samtools sort - -o ${name}
				"""
			};
			break;
		default: 
			error "Aligning with ${params.aligner} is not implemented. You can add it in 'align_reads' process"
			break;
	}
}

// DO not edit below
wasp_path = '/opt/WASP'

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

def filter_grouped_channel(ch) {
	ch.map(it -> tuple(it[0],
		it[1].findAll { f -> f[1] }.collect { element -> return element[0] })
		)
	filt_map
}

def get_container(file_name) {
  parent = file(file_name).parent
  container = "--bind ${parent}"
  if (file(file_name).exists()) {
	old_parent = file(file_name).toRealPath().parent
	if (old_parent != parent) {
		container += ",${old_parent}"
	}
  } 
  return container
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

process merge_snv_files {
	scratch true
	container "${params.container}"
	publishDir "${params.outdir}"

	input:
		path snv_files

	output:
		tuple path(name), path("${name}.tbi")

	script:
	name = "all_variants_stats.bed.gz"
	"""
	echo "${snv_files}" | tr ' ' '\n' | sort -t ":" -k1,1 -u > filelist.txt
	while read file; do
		zcat \$file | awk -v OFS='\t' -v name=`basename \$file | awk -F':' '{ print \$1 }'` '{print \$0,name}' >> variants.bed 
	done < filelist.txt

	sort-bed variants.bed | bgzip -c > ${name}
	tabix ${name}
	"""
}

process generate_h5_tables {
	scratch true
	publishDir "${params.outdir}/h5"
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

process split_reads {
	tag "${indiv_id}:${ag_number}"
	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)}"

	input:
		tuple val(ag_number), val(indiv_id), path(bam_file), path(bam_index_file), val(r_tag)

	output:
		tuple val(ag_number), val(indiv_id), val(r_tag), path(name), path("${name}.bai"), env(n_counts)

	script:
	name = "${ag_number}.${r_tag}.bam"
	pars = r_tag == 'pe' ? "-f 1" : "-F 1"
	"""
	samtools view -O bam -h ${pars} --reference ${params.genome_fasta_file} ${bam_file} > ${name}
	samtools index ${name}
	n_counts=\$(samtools view -c ${name})
	"""
}

process extract_to_remap_reads {
	tag "${ag_number}:${r_tag}"
	container "${params.container}"
	cpus 2
	scratch true

	input:
		tuple val(ag_number), val(indiv_id), val(r_tag), path(bam_file), path(bam_file_index), env(n_counts)
		path h5_tables

	output:
		tuple val(ag_number), val(r_tag), path(out_bam_file), path("${out_bam_file}.bai"), emit: bamfile
		tuple val(ag_number), val(r_tag), path(fasta1), path(fasta2), emit: fastq
	
	script:
	name = "${ag_number}.${r_tag}.rmdup"
	out_bam_file = "${name}.bam"
	fasta1 = "${name}.remap.fq1.gz"
	fasta2 = "${name}.remap.fq2.gz"
	if (r_tag == 'pe') {
		"""
		python3 ${wasp_path}/mapping/rmdup_pe.py \
			${bam_file} pe.reads.rmdup.bam

		samtools sort \
			-@${task.cpus} \
			-o ${out_bam_file} \
			-O bam \
			pe.reads.rmdup.bam

		samtools index ${out_bam_file}

		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
			--is_paired_end \
			--is_sorted \
			--output_dir ./ \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			${out_bam_file}
		"""
	} else {
		"""
		# an ugly hack to deal with repeated read names on legacy SOLEXA GA1 data
		python3 $moduleDir/bin/hash_se_reads.py ${bam_file} se.hashed.bam

		python3 ${wasp_path}/mapping/rmdup.py \
			se.hashed.bam  se.reads.rmdup.bam
		
		samtools sort \
			-@${task.cpus} \
			-o ${out_bam_file} \
			-O bam \
			se.reads.rmdup.bam
		
		samtools index ${out_bam_file}

		### Creates 3 following files:
		### se.reads.rmdup.sorted.to.remap.bam (reads to remap)
		### se.reads.rmdup.sorted.keep.bam (reads to keep)
		### se.reads.rmdup.sorted.remap.fq.gz (fastq file containing the reads with flipped alleles to remap)
		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
			--is_sorted \
			--output_dir ./ \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			${out_bam_file}
		
		mv ${name}.remap.fq.gz ${fasta1}
		ln -s ${fasta1} ${fasta2}
		"""
	}
}


process wasp_filter_reads {
	container "${params.container}"
	scratch true
	tag "${ag_number}:${r_tag}"

	input:
		tuple val(ag_number), val(r_tag), path(bam_file), path(initial_bam_file), path(initial_bam_file_index)
	
	output:
		tuple val(ag_number), path(name)

	script:
	name = "${ag_number}.${r_tag}.result.bam"
	"""
	python3 ${wasp_path}/mapping/filter_remapped_reads.py \
		${initial_bam_file} \
		${bam_file} \
		${name}
	"""
}

process merge_bam_files {
	container "${params.container}"
	scratch true
	tag "${ag_number}"

	input:
		tuple val(ag_number), path(bam_files)

	output:
		tuple val(ag_number), path(name), path("${name}.bai")

	script:
	name = "${ag_number}.merged.bam"
	if (bam_files.split(' ').size() >= 2)
		"""
		samtools merge -f reads.rmdup.original.bam \
			${bam_files}

		samtools sort \
			-@${task.cpus} \
			-o ${name} \
			reads.rmdup.original.bam
		samtools index ${name}
		"""
	else
		"""
		ln -s ${bam_files} ${name}
		samtools sort \
			-@${task.cpus} \
			-o ${name}  \
			reads.passing.bam
		samtools index ${name}
		"""
}

process calc_initial_read_counts {
	container "${params.container}"
	tag "${ag_number}"

	input:
		tuple val(ag_number), path(bam_file), path(bam_file_index)

	output:
		tuple val(ag_number), path(name), path("${name}.tbi")

	script:
	name = "${ag_number}.coverage.bed.gz"
	"""
	python3 $moduleDir/bin/count_tags_pileup.py ${filtered_sites_file} \
		${bam_file} \
		--only_coverage | bgzip -c > ${name}
	tabix ${name}
	"""
}

process count_reads {
	tag "${indiv_id}:${ag_number}"
	container "${params.container}"
	publishDir params.outdir + "/count_reads"

	input:
		tuple val(ag_number), val(indiv_id), path(filtered_sites_file), path(filtered_sites_file_index), path(bam_passing_file), path(bam_passing_file_index), path(rmdup_counts), path(rmdup_counts_index)

	output:
		tuple val(indiv_id), path(name), path("${name}.tbi")

	script:
	name = "${ag_number}.bed.gz"
	"""
	python3 $moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} ${bam_passing_file} \
		--original_dedup_cover ${rmdup_counts} \
		| sort-bed - | bgzip -c > ${name}
	tabix ${name}
	"""
}

process merge_by_indiv {
	publishDir "${params.outdir}/indiv_merged_files"
	tag "${indiv_id}"
	container "${params.container}"
	scratch true
	errorStrategy 'terminate'

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


workflow calcInitialReadCounts {
	take:
		data
	main:
		out = merge_bam_files(data) | calc_initial_read_counts
	emit:
		out
}

workflow waspRealigning {
	take:
		samples_aggregations
	main:
		h5_tables = generate_h5_tables().collect()

		sagr = samples_aggregations.map(it -> tuple(it[1], it[0], it[2], it[3]))
		
		indiv_ag_id_map = sagr.map(it -> tuple(it[0], it[1]))
		
		snps_sites = filter_variants(indiv_ag_id_map)
		merge_snv_files(snps_sites.map(it -> it[1]).collect(sort: true))
		
		samples = sagr.join(snps_sites, by: 0)
		r_tags = Channel.of('pe', 'se')
		split_rs = sagr
			| combine(r_tags)
			| split_reads
			| branch {
				files: it[5].toInteger() > 0
        		nodata: true
			}
		//tuple val(ag_number), val(indiv_id), val(r_tag), path(name), path("${name}.bai"), env(n_counts)
		// tuple val(ag_number), path(name)
		nodata = split_rs.nodata
			| map(it -> tuple(it[0], tuple(it[3], false)))

		to_remap_reads_and_initial_bam = extract_to_remap_reads(split_rs.files, h5_tables)

		dedup_bam = to_remap_reads_and_initial_bam.bamfile

		filtered_bam = to_remap_reads_and_initial_bam.fastq
			| align_reads
			| join(dedup_bam, by: [0, 1])
			| wasp_filter_reads
			| map(it -> tuple(it[0], tuple(it[1], true)))
		
		merged_out_bam = nodata	
			| mix(filtered_bam)
			| groupTuple(size: 2)
			| filter_grouped_channel
			| merge_bam_files

		initial_read_counts = nodata
			| mix(dedup_bam.map(it -> tuple(it[0], tuple(it[2], true))))
			| groupTuple(size: 2)
			| filter_grouped_channel
			| calcInitialReadCounts


		out = indiv_ag_id_map 
			| join(snps_sites)
			| join(merged_out_bam)
			| join(initial_read_counts)
			| count_reads
			| groupTuple()
			| merge_by_indiv
		
	emit:
		out
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
