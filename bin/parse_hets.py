import sys
import numpy as np

header = ["chrom", "start", "end", "variant_id", "dbsnp", "ref", "alt", "aa", "maf", "ard", "mu", "sigma", "n_ref", "n_alt", "n_total", "n_hets", "mean_rd"]
print('\t'.join(header))

for line in sys.stdin:
    fields = line.strip().split("\t")

    n = np.array(list(map(float, fields[8:-1:2])))
    ref = np.array(list(map(float, fields[9::2])))

    af = 'nan'
    if fields[7] != '.':
        af = sorted(map(float, filter(lambda x: x!='.', fields[7].split(','))), reverse=True)[1]
        af = '{:0.5f}'.format(af)
    fields[7] = af

    alt = n-ref

    ard = np.sum(ref)/np.sum(n)
    ard_mean = np.mean(ref/n)
    ard_sigma = np.std(ref/n)

    variant_id = f'{fields[0]}:{fields[2]}:{fields[4]}:{fields[5]}'
	
    if fields[3]=='.':
        fields[3] = f'{variant_id}\t{variant_id}'
    else:
        fields[3] = f'{variant_id}\t{fields[3]}'

    out = [
        ard, ard_mean, ard_sigma, np.sum(ref), np.sum(alt), np.sum(n), len(n), np.mean(n)
    ]

    print('\t'.join(fields[:8]) + "\t" + '\t'.join('{:0.4f}'.format(val) for val in out))
