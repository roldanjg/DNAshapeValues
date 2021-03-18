import pandas
import re
import os

fold = '/home/bgp01/tfbsshape/shapes/tomato_genome/'


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def get_files():
    for a, x, y in os.walk('.'):
        return y


def create_file(head, values, filename):
    with open(filename + '.wig', 'a+') as wiggle:
        wiggle.write('variableStep chrom={}'.format(head[1:]))
        for index, val in enumerate(values, 1):
            if not 'NA' in val:
                wiggle.write('{} {}\n'.format(index, val))


def split_chr(fasta):
    with open(fasta) as fa:
        headerchange = False
        for line in fa:
            if '>' in line:
                if headerchange:
                    create_file(header, numbers, fasta)
                    headerchange = False
                header = line
                numbers = []
            else:
                headerchange = True
                for i in line.split(','):
                    numbers.append(i if not '\n' in i else i[:-1])
        create_file(header, numbers, fasta)


with cd(fold):
    allfiles = get_files()
    for file in allfiles:
        if '.HelT' in file or '.MGW' in file or '.ProT' in file or '.Roll' in file:
            split_chr(file)
