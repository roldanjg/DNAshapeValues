import Bio.motifs as motifs
with open("meme51.xml") as handle:
    record = motifs.parse(handle, "meme")
motif = record[0]
motif_base = str(motif.counts)
with open('MA0566.1.pfm', 'a+') as transition:
     transition.write('>{}\t{}\n'.format('MA0566.1', 'MYC2'))
for line in motif_base.split('\n')[1:5]:
    base, values = line.strip().split(':')
    values = [str(round(float(i))) for i in values.split(' ') if i != '']
    with open('MA0566.1.pfm', 'a+') as transition:
        transition.write('{}  [ {} ]\n'.format(base, '  '.join(values)))