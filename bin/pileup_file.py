"""
Counts coverage of SNVs in the bam file
"""
from argparse import ArgumentParser
import pysam
import pandas as pd
import sys

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


def reads_to_dict(vars_file_path, bam_file_path, chrom):
    args = {} if chrom is None else {'reference': chrom} 
    with pysam.TabixFile(vars_file_path) as vars_file, pysam.AlignmentFile(bam_file_path, "rb") as sam_file: 
        for line in vars_file.fetch(*args):
            variant = SNV.from_str(line)
            if not variant.is_het:
                continue
            reads_1, reads_2, read_pairs = get_reads(variant, sam_file)
            yield variant, reads_1, reads_2, read_pairs


def main(var_file_path, bam_file_path, chrom):
    result = []
    for variant, _, _, read_pairs in reads_to_dict(var_file_path, bam_file_path, chrom):
        result.append([*variant.to_list(), len(read_pairs)])
    
    pd.DataFrame.from_records(result, columns=[*SNV.get_fields(), 'coverage']).to_csv(sys.stdout, 
    sep='\t', index=False, header=None)
        

if __name__ == '__main__':
    parser = ArgumentParser(description = "Count tags by allele")
    
    parser.add_argument("--chrom", dest = "chrom", type = str,
						default = None, help = "Use a specific contig/chromosome")
    parser.add_argument("var_file", metavar = "var_file", type = str,
						help = "Path to variant file (must have corresponding index)")

    parser.add_argument("remapped_bam_file", metavar = "remapped_bam_file", type = str, 
						help = "Path to BAM-format tag sequence file")

    args = parser.parse_args()

    main(args.var_file, args.remapped_bam_file, args.chrom)