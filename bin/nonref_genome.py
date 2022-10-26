#!/bin/env python

import sys
import shutil


import pyfaidx
import pysam


iupac = "XACMGRSVTWYHKDBN"
def get_iupac(ref, alts):
	i = iupac.find(ref)
	for alt in alts:
		j = iupac.find(alt)
		i = i | j
	return iupac[i]


orig_fasta_file=sys.argv[1]
out_fasta_file=sys.argv[2]
vcf_file=sys.argv[3]


print("Copying original fasta file...")

shutil.copy2(orig_fasta_file, out_fasta_file)
shutil.copy2(orig_fasta_file+".fai", out_fasta_file+".fai")

print('Changing nucleotides to iupac')
with pyfaidx.Fasta(out_fasta_file, mutable=True) as fasta:
    with pysam.VariantFile(vcf_file) as vcf:
        for rec in vcf.fetch():
            
            if any([len(i) > 1 for i in rec.alleles]):
                continue
            
            ref = rec.ref
            alts = rec.alts

            #if len(ref)>1 or len(alt)>1:
            #	continue

            ambig = get_iupac(ref, alts)
            try:
                fasta[rec.contig][rec.start] = ambig
            except KeyError as e:
                pass

                #print fasta[rec.contig][rec.start], ref, alt, ambig
