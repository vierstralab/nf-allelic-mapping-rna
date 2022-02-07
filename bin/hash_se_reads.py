#!/usr/bin/env python3

# adds a 4 character ASCII hash to read names; to deal with legacy SOLEXA GA1 data

import sys
import pysam
import random, string

original_file = sys.argv[1]
modified_file = sys.argv[2]

ofh = pysam.AlignmentFile(original_file)
nfh = pysam.AlignmentFile(modified_file, 'wb', template=ofh)

for read in ofh.fetch():

	if not read.is_paired:
		rand_tag = ''.join(random.choices(string.ascii_letters + string.digits, k=4))
		read.query_name += ':{}'.format(rand_tag)

	nfh.write(read)

ofh.close()
nfh.close()
