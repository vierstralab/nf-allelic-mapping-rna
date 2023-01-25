import pandas as pd
import sys

nucleotides = ('A', 'T', 'G', 'C')


def alt_str_has_single(alts_str):
    return sum([(len(alt) == 1) and (alt in nucleotides) for alt in alts_str.split(',')]) > 0


def main(snps, annotations):
    repeated = annotations[(annotations['topmed'] != '.') &
                           (annotations['topmed'].notna()) &
                           (annotations['alts'].apply(alt_str_has_single)) &
                           (annotations['ref'].isin(nucleotides))
    ].drop_duplicates(['#chr', 'end', 'ref', 'alts', 'topmed']).value_counts(['#chr', 'start', 'end', 'ref'])
    repeated = repeated[repeated > 1].reset_index()[['#chr', 'start', 'end', 'ref']]
    repeated = annotations.merge(repeated, on=['#chr', 'start', 'end', 'ref'])

    repeated = repeated.groupby(['#chr', 'start', 'end', 'ref'])['alts'].apply(
        lambda alts_on_same_pos:
        sum([alt_str_has_single(alts_str) for alts_str in alts_on_same_pos]) > 1
    )

    assert repeated.sum() == 0
    
    annotations = annotations[(annotations['topmed'] != '.') &
                           (annotations['topmed'].notna()) &
                           (annotations['alts'].apply(alt_str_has_single)) &
                           (annotations['ref'].isin(nucleotides))
    ]

    assert annotations.value_counts(['#chr', 'start', 'end', 'ref']).max() == 1

    merged = snps.merge(annotations, 
        on=['#chr', 'start', 'end', 'ref'],
        how='left')
    print(merged.columns)
    merged['RAF'] = merged['topmed'].apply(lambda x: '.' if pd.isna(x) or x == '.'
                                           else float(x.split(',')[0]))
    merged['AAF'] = merged.apply(
        lambda row:
        '.' if row['RAF'] == '.' else
        dict(zip(row['alts'].split(','), row['topmed'].split(',')[1:])).get(row['alt'], '.'),
        axis=1
    )
    return merged[['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id', 'AAF', 'RAF', 'FMR']]

if __name__ == '__main__':
    dbsnp_annotation = pd.read_table(sys.argv[1],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alts', 'topmed'])

    snps_to_annotate = pd.read_table(sys.argv[2], header=None, names=['#chr', 'start', 'end', 'ref', 'alts'])
    df = main(snps_to_annotate, dbsnp_annotation)
    df.to_csv(sys.argv[3], sep='\t', index=False, header=None)