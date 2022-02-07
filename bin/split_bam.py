#!/usr/bin/env python3

# Splits a BAM file into two files: one for single-end reads and one for pair-end reads
# it also adds a 4 character ASCII hash to single-end read names to deal with legacy SOLEXA GA1 data
# where there is no unique flowcell ID, which causes problems when merging multiple sequencing runs.

import sys
import pysam
import random, string

original_filename = sys.argv[1]
se_filename = sys.argv[2]
pe_filename = sys.argv[3]

original_bam = pysam.AlignmentFile(original_filename)
se_bam  = pysam.AlignmentFile(se_filename, 'wb', template=original_bam)
pe_bam = pysam.AlignmentFile(pe_filename, 'wb', template=original_bam)

for read in original_bam.fetch():

	if read.is_paired:
		pe_bam.write(read)
	else:
		rand_tag = ''.join(random.choices(string.ascii_letters + string.digits, k=4))
		read.query_name += ':{}'.format(rand_tag)
		se_bam.write(read)

original_bam.close()
se_bam.close()
pe_bam.close()

