import sys
import pandas as pd
import os


def get_snps_file_name(value, prefix):
    return os.path.join(prefix, value + '.snps.bed.gz')

def main(old_meta, output, file_prefix):
    df = pd.read_table(old_meta)
    df['snps_file'] = df['indiv_id'].apply(lambda x: get_snps_file_name(x, file_prefix))
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])

