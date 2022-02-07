#!/usr/bin/env nextflow

params.samples_file='/tmp/metadata.txt'
params.genotype_file='/net/seq/data/projects/regulotyping-h.CD3+/genotyping/output/calls/all.filtered.snps.annotated.vcf.gz'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'
params.outdir='output'

//DO NOT EDIT BELOW

nuclear_chroms = "$params.genome" + ".nuclear.txt"
genome_chrom_sizes_file="$params.genome"  + ".chrom_sizes"


Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.indiv_id, row.ag_number, row.bam_file, row.filtered_sites_file ) }
	.tap{ SAMPLES_AGGREGATIONS }

process generate_h5_tables {

	scratch true

	input:
		file vcf_file from file(params.genotype_file)
		file '*' from file("${params.genotype_file}.csi")
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

	snp2h5 --chrom chrom_sizes.txt.gz \
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
	module 'python/3.6.4'

	input:
	set val(indiv_id), val(ag_number), val(bam_file), val(filtered_sites_file) from SAMPLES_AGGREGATIONS

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
	set val(indiv_id), val(ag_number), val(filtered_sites_file), file("${ag_number}.bam"), file("${ag_number}.bam.bai"), file("${ag_number}.passing.bam"), file ("${ag_number}.passing.bam.bai") into REMAPPED_READS

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
		hash_se_reads.py se.bam se.hashed.bam

		## step 1 -- remove duplicates
		##
		python3 /home/jvierstra/.local/src/WASP/mapping/rmdup.py \
			se.hashed.bam  se.reads.rmdup.bam

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} \
			-o se.reads.rmdup.sorted.bam \
			-O bam \
			se.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		##
		python3 /home/jvierstra/.local/src/WASP/mapping/find_intersecting_snps.py \
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
		##

		python3 /home/solexa/stampipes/scripts/bwa/filter_reads.py \
			se.reads.remapped.bam \
			se.reads.remapped.marked.bam \
			${nuclear_chroms}

		## step 5 -- filter reads

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} -l0 se.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> se.reads.remapped.marked.filtered.bam

		python3 /home/jvierstra/.local/src/WASP/mapping/filter_remapped_reads.py \
			se.reads.rmdup.sorted.to.remap.bam \
			se.reads.remapped.marked.filtered.bam \
			se.reads.remapped.keep.bam

		merge_files="\${merge_files} se.reads.rmdup.sorted.bam"
		remapped_merge_files="\${remapped_merge_files} se.reads.remapped.keep.bam se.reads.rmdup.sorted.keep.bam"
	fi

	## paired-end
	if [[ \${n_pe} -gt 0 ]]; then
		
		## step 1 -- remove duplicates
		##
		python3 /home/jvierstra/.local/src/WASP/mapping/rmdup_pe.py \
			pe.bam pe.reads.rmdup.bam

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} \
			-o pe.reads.rmdup.sorted.bam \
			-O bam \
			pe.reads.rmdup.bam

		## step 2 -- get reads overlapping a SNV
		##
		python3 /home/jvierstra/.local/src/WASP/mapping/find_intersecting_snps.py \
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

		python3 /home/solexa/stampipes/scripts/bwa/filter_reads.py \
			pe.reads.remapped.bam \
			pe.reads.remapped.marked.bam \
			${nuclear_chroms}

		## step 5 -- filter reads

		samtools sort \
			-m ${task.memory.toMega()}M \
			-@${task.cpus} -l0 pe.reads.remapped.marked.bam \
		| samtools view -b -F 512 - \
		> pe.reads.remapped.marked.filtered.bam

		python3 /home/jvierstra/.local/src/WASP/mapping/filter_remapped_reads.py \
			pe.reads.rmdup.sorted.to.remap.bam \
			pe.reads.remapped.marked.filtered.bam \
			pe.reads.remapped.keep.bam

		merge_files="\${merge_files} pe.reads.rmdup.sorted.bam"
		remapped_merge_files="\${remapped_merge_files} pe.reads.remapped.keep.bam pe.reads.rmdup.sorted.keep.bam"
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

	samtools merge -f \
		reads.passing.bam \
		\${remapped_merge_files}

	samtools sort \
		-m ${task.memory.toMega()}M \
		-@${task.cpus} \
		-o reads.passing.sorted.bam  \
		reads.passing.bam 

	# todo: merge dedupped se and pe reads
	mv reads.rmdup.sorted.bam ${ag_number}.bam
	samtools index ${ag_number}.bam

	mv reads.passing.sorted.bam ${ag_number}.passing.bam
	samtools index ${ag_number}.passing.bam

	"""
}

process count_reads {
	tag "${indiv_id}:${ag_number}"

	publishDir params.outdir + "/count_reads", mode: 'symlink'

	module 'python/3.6.4'

	input:
	set val(indiv_id), val(ag_number), val(filtered_sites_file), file(bam_all_file), file(bam_all_index_file), file(bam_passing_file), file(bam_passing_index_file) from REMAPPED_READS

	output:
	set val(indiv_id), val(ag_number), file("${ag_number}.bed.gz"), file("${ag_number}.bed.gz.tbi") into COUNT_READS_LIST, COUNT_READS_FILES

	script:
	"""
	count_tags_pileup.py \
		${filtered_sites_file} ${bam_all_file} ${bam_passing_file} \
	| sort-bed - | bgzip -c > ${ag_number}.bed.gz

	tabix -p bed ${ag_number}.bed.gz
	"""
}

// ag num, indiv, read file
COUNT_READS_LIST
	.map{ [it[1], it[0], it[2].name].join("\t") }
	.collectFile(
		name: 'sample_map.tsv', 
		newLine: true
	)
	.first()
	.set{RECODE_VCF_SAMPLE_MAP_FILE}

COUNT_READS_FILES
	.flatMap{ [file(it[2]), file(it[3])] }
	.set{ COUNT_READS_ALL_FILES }

process recode_vcf {
	publishDir params.outdir + "/vcf", mode: 'symlink'

	module 'python/3.6.4'

	input:
	file 'sample_map.tsv' from RECODE_VCF_SAMPLE_MAP_FILE
	file '*' from COUNT_READS_ALL_FILES.collect()

	file vcf_file from file("${params.genotype_file}")
	file '*' from file("${params.genotype_file}.csi")

	output:
	file 'allele_counts.vcf.gz*'

	script:
	"""
	recode_vcf.py \
		${vcf_file} \
		sample_map.tsv \
		allele_counts.vcf

	bgzip -c allele_counts.vcf > allele_counts.vcf.gz
	bcftools index allele_counts.vcf.gz

	"""
}
