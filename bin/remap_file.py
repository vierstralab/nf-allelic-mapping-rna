import pysam
import sys
import os
import pandas as pd
header = ['chr', 'start', 'end', 'ID', 'ref', 'alt', 'maf',
 'gt', 'ref_reads', 'alt_reads', 'original', 'failed_mapping',
 'n_failed_genotyping', 'n_failed_bias']

df_new = pd.read_table(sys.argv[1], header=None, names=header)
df_old = pd.read_table(sys.argv[2], header=None, names=header)

assert len(df_new.index) == len(df_old.index)
result_df = df_new.merge(df_old, on=['chr', 'start', 'end', 'ID', 'ref', 'alt', 'maf',
 'gt'])
assert len(result_df.index) == len(df_new.index)

result_df[['ref_reads', 'alt_reads', 'n_failed_genotyping', 'n_failed_bias']] = \
    result_df[['ref_reads_x', 'alt_reads_x', 'n_failed_genotyping_x', 'n_failed_bias_x']]

result_df['original'] = result_df['original_y'] - result_df['failed_mapping_y']
if (result_df['original'] < 0).sum() != 0:
    print(result_df[result_df['original'] <= 0])
    raise AssertionError

result_df['failed_mapping'] = result_df['failed_mapping_x'] - result_df['failed_mapping_y']
assert (result_df['failed_mapping'] < 0).sum() == 0

result_df[header].to_csv(sys.stdout, sep='\t', header=None, index=False)