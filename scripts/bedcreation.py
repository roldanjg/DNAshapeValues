import pandas
import re
import sys
fasta = sys.argv[1]
finalname = str(sys.argv[2])
def correction(line):

    regular = re.match(r"(.+):([0-9]+)-([0-9]+)", line)
    return regular.group(1), regular.group(2), regular.group(3)


def create_files(filetype, dataframe):
    with open(filetype+'.bed', 'w') as finalbed:
        for index, row in dataframe.iterrows():
            finalbed.write('{}\t{}\t{}\t{}\n'.format(row.chr, row.start, row.end, row.id))



with open(fasta, 'r') as totalfastafile:
    df = pandas.read_csv(totalfastafile, sep='>', names=['sequence', 'id'])
dfid = df['id'].dropna().reset_index(drop=True)
dfsequence = df['sequence'].dropna().reset_index(drop=True)
dfasta = pandas.merge(dfid, dfsequence, left_index=True, right_index=True)
df = dfasta
df['chr'], df['start'], df['end'] = zip(*df['id'].map(correction))

create_files(finalname, df)