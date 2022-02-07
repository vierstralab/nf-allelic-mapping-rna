import sys
import numpy as np

header = ["variant_id", "dbsnp", "ref", "alt",  "indiv_id", "celltype", "nref", "nalt"]
print('\t'.join(header))

def unpack(s):
    sample_id, ard = s.split(':')
    total, ref = ard.split('/')
    total = int(total)
    ref = int(ref)

    split_index = sample_id.rfind('_')
    indiv = sample_id[:split_index]
    cell_type = sample_id[split_index+1:]

    return indiv, cell_type, total, ref

for line in sys.stdin:
    fields = line.strip().split("\t")

    af = 'nan'
    if fields[7] != '.':
        af = sorted(map(float, filter(lambda x: x!='.', fields[7].split(','))), reverse=True)[1]
        #af = float(fields[7].split(',')[0])
        af = '{:0.5f}'.format(af)
    fields[7] = af


    variant_id = f'{fields[0]}:{fields[2]}:{fields[4]}:{fields[5]}'
	
    if fields[3]=='.':
        fields[3] = f'{variant_id}\t{variant_id}'
    else:
        fields[3] = f'{variant_id}\t{fields[3]}'


    samples = fields[8:]
    
    ratio_dict = {}

    for sample in samples:
        indiv, cell_type, total, ref = unpack(sample)

        ratio = np.log((ref+1) / (total-ref+1))
        #ratio = ref / total

        if cell_type not in ratio_dict:
            ratio_dict[cell_type] = list()
    
        ratio_dict[cell_type].append(f'{ratio:0.3f}')


    out = []
    for k, v in ratio_dict.items():
        out.append(f'{k}:{",".join(map(str, v))}')

    print('\t'.join(fields[:8]) + '\t' + ';'.join(out))