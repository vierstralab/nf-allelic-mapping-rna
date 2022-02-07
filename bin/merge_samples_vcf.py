#!/usr/bin/env python3

import sys
import logging
from argparse import ArgumentParser

from collections import namedtuple

import pysam
import numpy as np

import pandas as pd

def parse_options(args):

    parser = ArgumentParser(description = "Combined allelic read depths per sample into a larger VCF file")

    parser.add_argument("variant_file", metavar = "variant_file", type = str,
                        help = "Path to VCF-format read depth file.")

    parser.add_argument("sample_map_file", metavar = "sample_map_file", type = str,
                        help = "Sample to individual mapping file")

    parser.add_argument("outfile_per_group", metavar = "outfile_per_group", type = str,
                        help = "Output VCF file")

    parser.add_argument("--chrom", metavar = "chrom", type = str,
                    default=None, help = "Restrict to a specific chromosome")

    return parser.parse_args(args)


# Make a new VCF file
def make_vcf_header_from_template(template_header, samples):
    header=pysam.VariantHeader()
    for record in template_header.records:
        header.add_record(record)

    # hack because the clearing the formats doesn't actually delete them "behind the scenes"
    header.formats.clear_header()
    header=header.copy()

    header.formats.add('GT', 1, "String", "Genotype")
    header.formats.add('AD', 2, "Integer", "Allele read depth")
    header.formats.add('RD', 1, "Integer", "Read depth (unfiltered)")
    header.formats.add('ARD', 1, "Integer", "Read depth (passing filters)")
    header.formats.add('FM', 1, "Integer", "Failed  - mapping")
    header.formats.add('FG', 1, "Integer", "Failed - discordant genotypes")
    header.formats.add('FB', 1, "Integer", "Failed - 5' proximity bias")
    header.formats.add('FMR', 1, "Float", "Failed mapping rate")

    for sample in samples:
        header.add_sample(sample)

    return header

class variant_data:
    def __init__(self, gt=(None, None), ad=[0,0], ard=0, rd=0, fm=0, fg=0, fb=0, phased=False, **kwargs):
        self.gt = gt
        self.ad = list(ad)
        self.ard = ard
        self.rd = rd
        self.fm = fm
        self.fg = fg
        self.fb = fb

        self.fmr = 0 if self.rd == 0 else self.fm/self.rd
        self.phased = phased

    def __add__(self, other):

        if self.gt != other.gt:
            print(self.gt, other.gt)
            raise ValueError("Genotypes do not match!")

        self.ad[0] += other.ad[0]
        self.ad[1] += other.ad[1]
        self.ard += other.ard
        self.rd += other.rd
        self.fm += other.fm
        self.fg += other.fg
        self.fb += other.fb
        
        self.fmr = 0 if self.rd == 0 else self.fm/self.rd
        self.phased = self.phased & other.phased

        return self

    def _asdict(self):
        return {
            'GT': self.gt,
            'AD': self.ad,
            'ARD': self.ard,
            'RD': self.rd,
            'FM': self.fm,
            'FG': self.fg,
            'FB': self.fb,
            'FMR': self.fmr,
            'phased': self.phased,
        }
    
    @classmethod
    def _fromdict(cls, x):
        return cls(**{k.lower(): v for k, v in x.items()})

def main(argv = sys.argv[1:]):

    args = parse_options(argv)

    # samples=[]
    # samples_group_id={}

    # with open(args.sample_map_file) as f:
    #     for line in f:
    #         (genotype_sample_id, sample_id, sample_cell_type) = line.strip().split("\t")

    #         samples.append(sample_id)
    #         samples_group_id[sample_id] = f'{genotype_sample_id}'
    #         #samples_group_id[sample_id] = f'{genotype_sample_id}_{cell_type}'

    samples = pd.read_table(args.sample_map_file, header=0, dtype={'ag_number': str})
    samples.set_index('ag_number', inplace=True)

    samples['group_id'] = samples['indiv_id'] + '_' + samples['cell_type']

    samples_group_id = samples['group_id'].to_dict()

    group_ids = list(set(samples_group_id.values()))

    infile = pysam.VariantFile(args.variant_file, mode='r', ignore_truncation=True)

    outfile_per_group = pysam.VariantFile(args.outfile_per_group, mode='w', 
        header = make_vcf_header_from_template(infile.header, group_ids))
    outfile_per_group.header.add_meta('ARD_command', 'merge_samples_vcf.py ' + ' '.join(['%s=%s' %(k, v) for k, v in vars(args).items()]))

    for var in infile.fetch(contig=args.chrom):
        vd_grouped = {}
  
        for sample in samples.index:

            group_id = samples_group_id[sample]
            if group_id in vd_grouped:
                # Raises 'ValueError' if genotypes don't match!
                vd_grouped[group_id] += variant_data._fromdict(var.samples[sample])
            else:
                vd_grouped[group_id] = variant_data._fromdict(var.samples[sample])

        samples_format = [vd_grouped.get(group_id, variant_data())._asdict() for group_id in group_ids]
 
        outvar = outfile_per_group.new_record(contig=var.contig, start=var.start, stop=var.stop, alleles=var.alleles, id=var.id, qual=var.qual, filter=var.filter, info=var.info, samples=samples_format)
        outfile_per_group.write(outvar)

    outfile_per_group.close()
    infile.close()

    return 0
    
if __name__ == "__main__":
    sys.exit(main())

