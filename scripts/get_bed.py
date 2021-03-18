import os
import pandas
import re
fasta = 'MpoMYCH1_51bp.fa'
bed = 'MpoMYCH1_51bp.bed'

with cd('C:/Users/JOAQUINGR/PycharmProjects/untitled/mldata/marchantia/MYC2/'):
    with open(fasta, 'r') as totalfastafile:
        df = pandas.read_csv(totalfastafile, sep='>', names=['sequence', 'ide'])
    dfid = df['ide'].dropna().reset_index(drop=True)
    dfsequence = df['sequence'].dropna().reset_index(drop=True)
    dfasta = pandas.merge(dfid, dfsequence, left_index=True, right_index=True)

    with open(bed, 'r') as totalbedfile:
            dfbed = pandas.read_csv(totalbedfile, sep='\t', names=['chr', 'start', 'end', 'id'])
    df = pandas.concat([dfbed, dfasta], axis=1)
    del dfid, dfsequence, dfasta, dfbed
    with open(bed+ '.bed', 'w') as finalbed:
        with open(fasta+ '.fa', 'w') as finalfa:
            for index, row in df.iterrows():
                finalbed.write('{}\t{}\t{}\t{}\n'.format(row.chr, row.start, row.end, row.ide))
                finalfa.write('>{}\n{}\n'.format(row.ide, row.sequence))