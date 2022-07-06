import pysam
import sys
import os

from pileup_file import SNV


def main(in_file):
    with pysam.TabixFile(in_file) as vars_file:
        for line in vars_file.fetch():
            split_line = line.strip('\n').split('\t')
            variant = SNV(split_line[:-6])
            n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias = split_line[-6:]
            print(line, n_ref, n_alt)
            assert n_original_reads == n_alt + n_ref + n_failed_bias + n_failed_genotyping + n_failed_mapping
            if min(n_ref, n_alt) < 5:
                continue
            if variant.maf < 0.05:
                continue
            print('\t'.join(map(str, [*split_line[:6], n_ref, n_alt])))


if __name__ == '__main__':
    main(sys.argv[1])