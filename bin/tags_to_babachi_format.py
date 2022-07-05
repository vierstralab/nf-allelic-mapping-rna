import pandas as pd
import sys
import os

from pileup_file import SNV


def main(in_file):
    print(os.path.exists(in_file))
    df = pd.read_table(in_file)
    columns = SNV.get_fields()[:6] + ['ref_counts', 'alt_counts']
    df[columns].to_csv(sys.stdout, sep='\t', header=None, index=False)


if __name__ == '__main__':
    main(sys.argv[1])