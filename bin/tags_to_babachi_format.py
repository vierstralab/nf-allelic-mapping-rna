import pysam
import sys
import os

from count_tags_pileup import SNV


def main(in_file):
    name = os.path.basename(in_file).split('.')[0]
    with pysam.TabixFile(in_file) as vars_file:
        for line in vars_file.fetch():
            split_line = line.strip('\n').split('\t')
            variant = SNV(split_line[:-6])
            n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias = map(int, split_line[-6:])
            assert n_original_reads == n_alt + n_ref + n_failed_bias + n_failed_genotyping + n_failed_mapping
            print('\t'.join(map(str, [*split_line[:6], n_ref, n_alt, name, variant.aaf, variant.raf, n_failed_mapping / n_original_reads])))


if __name__ == '__main__':
    main(sys.argv[1])
