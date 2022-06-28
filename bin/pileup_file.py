"""
Counts coverage of SNVs in the bam file
"""

import sys
import logging

from argparse import ArgumentParser

import pysam

class SNV:
    """chrom, start, end, id, ref, alt, gt, extra fields
        GT encoded as either 0/1 or with pipe 0|0
    """
    
    __class_fields = ['contig', 'start', 'end', 'id', 'ref', 'alt', 'gt']
    def __init__(self, fields):
        for field_name, field_value in zip(self.__class_fields, fields):
            setattr(self, field_name, field_value)
        self.start = int(self.start)
        self.end = int(self.end)
        self.is_het = sum(map(int, self.gt.replace('|','/').split('/'))) == 1
       
    def to_list(self):
        return [getattr(self, field) for field in self.__class_fields]

    def __repr__(self):
        return '\t'.join(map(str, self.to_list()))

    @staticmethod
    def from_str(line):
        return SNV(line.strip('\n').split('\t'))
    
    @classmethod
    def get_fields(cls):
        return cls.__class_fields

def get_reads(variant, sam_file):

	reads_1 = {}
	reads_2 = {}

	# Go into BAM file and get the reads
	for pileupcolumn  in sam_file.pileup(variant.contig, variant.start, variant.start+1, maxdepth=10000, truncate=True, stepper="nofilter"):

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
    with pysam.TabixFile(vars_file_path) as vars_file, pysam.AlignmentFile(bam_file_path, "rb") as sam_file: 
        for line in vars_file.fetch(reference=chrom):
            variant = SNV.from_str(line)
            if not variant.is_het:
                continue
            reads_1, reads_2, read_pairs = get_reads(variant, sam_file)
            yield variant, reads_1, reads_2, read_pairs


def main(var_file_path, bam_file_path, chrom):
    print('\t'.join(SNV.get_fields()))
    for variant, _, _, read_pairs in reads_to_dict(var_file_path, bam_file_path, chrom):
        print(str(variant), len(read_pairs), sep='\t')
        

if __name__ == '__main__':
    main()