#!/usr/bin/env python3
# Jeff Vierstra 2018
# TODO:
# --add filters/etc. as option
import sys
import logging

from argparse import ArgumentParser
import pysam

logging.basicConfig(stream = sys.stderr, level='WARNING')


class SNV:
    """chrom, start, end, id, ref, alt, maf, gt
        GT encoded as either 0/1 or with pipe 0|0
    """
    
    __class_fields = ['contig', 'start', 'end', 'id', 'ref', 'alt', 'maf', 'gt', 'gq', 'n_original_reads']
    def __init__(self, fields):
        for field_name, field_value in zip(self.__class_fields, fields):
            setattr(self, field_name, field_value)
        self.start = int(self.start)
        self.end = int(self.end)
        self.maf = float(self.maf)
        self.is_het = sum(map(int, self.gt.replace('|','/').split('/'))) == 1
        self.n_original_reads = int(self.n_original_reads)
       
    def to_list(self):
        return [getattr(self, field) for field in self.__class_fields]

    def __repr__(self):
        return '\t'.join(map(str, self.to_list()))

    @classmethod
    def from_str(cls, line: str):
        return cls(line.strip('\n').split('\t'))
    
    @classmethod
    def get_fields(cls):
        return cls.__class_fields

def reads_to_dict(vars_file_path, bam_file_path, chrom):
    args = {} if chrom is None else {'reference': chrom} 
    with pysam.TabixFile(vars_file_path) as vars_file, pysam.AlignmentFile(bam_file_path, "rb") as sam_file: 
        for line in vars_file.fetch(*args):
            variant = SNV.from_str(line)
            if not variant.is_het:
                continue
            reads_1, reads_2, read_pairs = get_reads(variant, sam_file)
            yield variant, reads_1, reads_2, read_pairs


def get_reads(variant, sam_file):

	reads_1 = {}
	reads_2 = {}

	# Go into BAM file and get the reads
	for pileupcolumn  in sam_file.pileup(variant.contig, variant.start, variant.end,
     maxdepth=10000, truncate=True, stepper="nofilter"):

		for pileupread in pileupcolumn.pileups:

			if pileupread.is_del or pileupread.is_refskip:
				continue

			if pileupread.alignment.is_read1:
				reads_1[pileupread.alignment.query_name] = pileupread
			else:
				reads_2[pileupread.alignment.query_name] = pileupread

	# All reads that overlap SNP; unqiue set
	read_pairs = set(reads_1.keys()) | set(reads_2.keys())

	return reads_1, reads_2, read_pairs


def parse_options(args):

	parser = ArgumentParser(description = "Count tags by allele")

	parser.add_argument("--chrom", dest = "chrom", type = str,
						default = None, help = "Use a specific contig/chromosome")

	parser.add_argument("var_file", metavar = "var_file", type = str,
						help = "Path to variant file (must have corresponding index)")

	parser.add_argument("remapped_bam_file", metavar = "remapped_bam_file", type = str, 
						help = "Path to BAM-format tag sequence file")

	return parser.parse_args(args)

class GenotypeError(Exception):
	pass

class DiploidError(Exception):
	pass

class ReadBiasError(Exception):
	pass

class ReadAlignmentError(Exception):
	pass

class ReadGenotypeError(Exception):
	pass

def get_5p_offset(pileupread):
	"""
	Returns position of variant relative to 5' of read
	"""
	if pileupread.query_position is None: # pileup overlaps deletion 
		return None
	elif pileupread.alignment.is_reverse:
		return pileupread.alignment.query_length-pileupread.query_position
	else:
		return pileupread.query_position+1

def get_base_quality(pileupread):
	"""
	Returns base call quality at variant position
	"""
	return pileupread.alignment.query_qualities[pileupread.query_position]
					

def check_bias(pileupread, offset=3, baseq=20):

	if pileupread is None:
		return True

	if get_5p_offset(pileupread)<=offset:
		raise ReadBiasError()

	# if get_base_quality(pileupread)<baseq:
	# 	raise ReadBiasError()

	return True

def get_base_call(pileupread):

	if pileupread is None:
		return None

	if pileupread.query_position is None:
		return None
	else:
		return pileupread.alignment.query_sequence[pileupread.query_position]

def check_alleles(pileupread, ref_allele, nonref_allele):

	if pileupread is None:
		return True

	# if pileupread.alignment.mapping_quality<30:
	# 	raise ReadAlignmentError()
	
	read_allele = get_base_call(pileupread)
	if read_allele != ref_allele and read_allele != nonref_allele:
		return ReadGenotypeError()

	# if read_allele == ref_allele:
	# 	num_permitted_mismatches = 1 
	# elif read_allele == nonref_allele:
	# 	num_permitted_mismatches = 2 
	# else:
	# 	return ReadGenotypeError()

	# mismatches = int(pileupread.alignment.get_tag("XM", with_value_type=False))
	# if mismatches > num_permitted_mismatches:
	# 	raise ReadAlignmentError()

	# if re.search("[^ACGT]", pileupread.alignment.query_sequence):
	# 	raise ReadAlignmentError()
	# 	# raise AlignmentError("Ambiguous base calls within read (not matching {A, C, G, T})")

	# if re.search("[HSPDI]", pileupread.alignment.cigarstring):
	# 	raise ReadAlignmentError()
	# 	# raise AlignmentError("Deletions/indels within read")

	return True

def check_reads(reads_1, reads_2, unique_reads, ref, alt):
	n_ref = n_alt = n_failed_bias = n_failed_genotyping = 0
	for read in unique_reads:
		try:

			read1 = reads_1.get(read, None)
			check_alleles(read1, ref, alt) # only returns true if read exists
			check_bias(read1) # only returns true if read exists
			read1_allele = get_base_call(read1) # returns None if read doesn't exist
			
			read2 = reads_2.get(read, None)
			check_alleles(read2, ref, alt) # only returns true if read exists
			check_bias(read2) # only returns true if read exists
			read2_allele = get_base_call(read2) # returns None if read doesn't exist

			read_allele = read1_allele or read2_allele

			# No ba errors
			if read_allele == ref:
				n_ref += 1
			elif read_allele == alt:
				n_alt += 1
			else:
				raise ReadGenotypeError()

		except ReadBiasError as e:
			n_failed_bias += 1
			logging.debug("Failed bias: " + read)
			continue
		except ReadGenotypeError as e:
			n_failed_genotyping += 1
			logging.debug("Failed genotyping: " + read)
			continue
	return n_ref, n_alt, n_failed_bias, n_failed_genotyping

	
def main(argv = sys.argv[1:]):

	args = parse_options(argv)
	for variant, reads_1, reads_2, read_pairs in reads_to_dict(args.var_file, args.remapped_bam_file, args.chrom):
		n_remapped_reads = len(read_pairs)
		n_ref, n_alt, n_failed_bias, n_failed_genotyping = check_reads(reads_1, reads_2,
															read_pairs, variant.ref, variant.alt)
		variant_str = str(variant)
		n_original_reads = variant.n_original_reads
		
		n_failed_mapping = n_original_reads - n_remapped_reads

		print(variant_str, n_ref, n_alt, n_original_reads, n_failed_mapping,
		 n_failed_genotyping, n_failed_bias, sep='\t')
    
if __name__ == "__main__":
    main()

