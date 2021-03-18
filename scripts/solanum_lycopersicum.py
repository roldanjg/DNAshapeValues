from Bio import SeqIO
from Bio.SeqUtils import GC
import Bio.motifs as motifs
import statistics
import urllib.request
import gzip
import shutil
import csv
import subprocess
import re
import os
from pathlib import Path
import time
import pandas
import numpy as np
from sklearn.model_selection import StratifiedKFold


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class generate_training_files:

    def __init__(self, species, tf, peaksize, truefa, truebed, memexml, memetxt, shape, url, script):
        self.species_name = species
        self.tf_name = tf
        self.peaksize = peaksize
        self.true_peaks = truefa
        self.true_bed = truebed
        self.memexml = memexml
        self.memetxt = memetxt
        self.shape = shape
        self.url = url
        self.script_base = script
        self.pfm = species + tf + '_' + peaksize + '.pfm'
        self.false_peaks_biasaway = self.true_peaks + 'background.biasaway.fa'
        self.randomsequences = 'random.' + species + '.bed'
        self.randomfasta = 'random.' + species + '.fa'
        self.genome = None
        self.n_samples = None
        self.gc_total = None
        self.fistorder = None
        self.detailed = None

    def create_pfm_from_meme(self, memexml, nameforpfm, tf_name, species_name):
        with open(memexml) as handle:
            record = motifs.parse(handle, "meme")
        motif = record[0]
        motif_base = str(motif.counts)
        with open(nameforpfm, 'w') as transition:
            transition.write('>{}\t{}\n'.format(tf_name, species_name))
            for line in motif_base.split('\n')[1:5]:
                base, values = line.strip().split(':')
                values = [str(round(float(i))) for i in values.split(' ') if i != '']
                transition.write('{}  [ {} ]\n'.format(base, '  '.join(values)))

    def create_folders(self):

        base = '/home/bgp01/tfbsshape/data/'

        shapes_vals = os.path.join(base, self.species_name, "shape_vals")
        Path(shapes_vals).mkdir(parents=True, exist_ok=True)

        Path(os.path.join(base, self.species_name, self.tf_name, self.peaksize)).mkdir(parents=True, exist_ok=True)

        self.learning = os.path.join(base, self.species_name, self.tf_name, self.peaksize, "learning")
        Path(self.learning).mkdir(parents=True, exist_ok=True)

        self.pathshuffle = os.path.join(self.learning, 'shuffle')
        Path(self.pathshuffle).mkdir(parents=True, exist_ok=True)

        self.pathcg = os.path.join(self.learning, 'cg')
        Path(self.pathcg).mkdir(parents=True, exist_ok=True)

        self.common_files = os.path.join(base, self.species_name, self.tf_name, self.peaksize, "common_files")
        Path(self.common_files).mkdir(parents=True, exist_ok=True)

        #         this is a folder where you might put all the data that initially is going to be clasified (init args)
        base = '/home/bgp01/tfbsshape/data/clasification'
        # move fasta file

        path1 = os.path.join(base, self.true_peaks)
        path2 = os.path.join(self.common_files, self.true_peaks)
        os.rename(path1, path2)

        # move bed file

        path1 = os.path.join(base, self.true_bed)
        path2 = os.path.join(self.common_files, self.true_bed)
        os.rename(path1, path2)

        # move pfm file

        path1 = os.path.join(base, self.pfm)
        path2 = os.path.join(self.common_files, self.pfm)
        os.rename(path1, path2)

        # move self.fistorder file

        path1 = os.path.join(base, self.memexml)
        path2 = os.path.join(self.common_files, self.memexml)
        os.rename(path1, path2)

        # move detailed file

        path1 = os.path.join(base, self.memetxt)
        path2 = os.path.join(self.common_files, self.memetxt)
        os.rename(path1, path2)

        # move shapevals file
        try:
            path1 = os.path.join(base, self.shape)
            path2 = os.path.join(shapes_vals, self.shape)
            os.rename(path1, path2)
        except:
            print('''data for shape values is already created, or is not in forlder, 
                    if this is the first time you create this specie, clt+C and check it''')

    def download_genome(self, url):

        file_name = url.split('/')[-1]
        self.genome = file_name[:-3]
        print(file_name)
        print(self.genome)
        urllib.request.urlretrieve(url, file_name)
        if '.gz' in file_name:
            with gzip.open(file_name, 'rb') as f_in:
                with open(self.genome, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    def calcularte_CG_content(self, peaks_fa_file):
        gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse(peaks_fa_file, "fasta"))
        self.gc_total = statistics.median(gc_values)

    def generate_false_peaks_cg(self):

        subprocess.call(['/home/bgp01/webproyect/programs/samtools-1.10/samtools',
                         'faidx', self.genome])
        with open(self.genome + '.fai') as index:
            for line in index:
                chr, number, *args = line.split('\t')
                with open('index.genome', 'a+') as genome:
                    genome.write('{}\t{}\n'.format(chr, number))
        f = open(self.randomsequences, "w")
        subprocess.call(['/home/bgp01/webproyect/programs/bedtools2/bedtools',
                         'random', '-l', '51', '-g', 'index.genome', '-n', '1300000'], stdout=f)
        subprocess.call(['/home/bgp01/webproyect/programs/bedtools2/bedtools',
                         'getfasta', '-name', '-s', '-fi', self.genome, '-bed', self.randomsequences,
                         '-fo', self.randomfasta])

    def generate_false_peaks_biasaway(self):
        #         subprocess.call(['tmux', 'new-session', '-d', '-s', speciename+'shape'])
        #         subprocess.call(['tmux', 'send-keys', '-t', speciename+'shape' + ':0', 'conda ', 'activate ',
        #                          'tfbsshape3', 'C-m'])
        f = open(self.false_peaks_biasaway, "w")
        subprocess.call(['biasaway', 'd', '-f', self.true_peaks], stdout=f)
        print('done')

    def split_files_cv(self, fasta, top_peaks_pwm, num_splits, name, bed=False):

        def correction(line):
            try:
                regular = re.match(r'([0-9]+)::(.+):([0-9]+)-([0-9]+)\([+,-]\)', line)
                return regular.group(2), regular.group(3), regular.group(4)
            except:
                print('this line dont match the regular expresion for bed cg%: ', line)

        def create_files(filetype, n_CV, dataframe):
            with open(filetype + name + str(n_CV) + '.bed', 'w') as finalbed:
                with open(filetype + name + str(n_CV) + '.fa', 'w') as finalfa:
                    for index, row in dataframe.iterrows():
                        finalbed.write('{}\t{}\t{}\t{}\n'.format(row.chr, row.start, row.end, row.id))
                        finalfa.write('>{}\n{}\n'.format(row.id, row.sequence))

        with open(fasta, 'r') as totalfastafile:
            df = pandas.read_csv(totalfastafile, sep='>', names=['sequence', 'id'])
        dfid = df['id'].dropna().reset_index(drop=True)
        dfsequence = df['sequence'].dropna().reset_index(drop=True)
        dfasta = pandas.merge(dfid, dfsequence, left_index=True, right_index=True)

        if 'foreground' in name:

            with open(bed, 'r') as totalbedfile:
                dfbed = pandas.read_csv(totalbedfile, sep='\t', names=['chr', 'start', 'end', 'id'])
            df = pandas.merge(dfbed, dfasta, on="id")
            del dfid, dfsequence, dfasta, dfbed
            self.n_samples, _ = df.shape


        elif 'cg' in name and self.n_samples != None:

            if self.n_samples == None:
                raise (print('total number of peaks in foreground is not defined, run foreground before background'))
            df = dfasta
            # create all the bed collumns based on >names from random sequences fasta in the dataframe
            df['chr'], df['start'], df['end'] = zip(*df['id'].map(correction))
            true_mean_gc = round(self.gc_total)
            print(true_mean_gc, 'rounded total mean')
            # extract from the 1000000 list only the ones that has same GC content as peaks mean
            df = df.loc[df['sequence'].apply(lambda x: round(GC(x))) == true_mean_gc].reset_index(drop=True)
            print(df.shape, 'hola')
            print(df.head())
            # keep only as much sequences as in seq peaks file
            df = df.loc[:self.n_samples - 1]
            print(df.shape, 'qtal')
            print(df.head())

        df['class'] = np.hstack(([1] * self.n_samples))
        # delete the peaks used to create the PWM
        df = df.loc[int(top_peaks_pwm - 1):].reset_index(drop=True)
        print(df.shape)
        skf = StratifiedKFold(n_splits=num_splits, shuffle=True)
        skf.get_n_splits(df, df['class'])
        for count, (train_index, test_index) in enumerate(skf.split(df, df['class'])):
            X_train, X_test = df.loc[train_index], df.loc[test_index]
            create_files('training_', count, X_train)
            create_files('testing_', count, X_test)

    def data_clasification(self):
        back = "background"
        fore = 'foreground'
        all_files = []

        for base in [self.pathshuffle, self.pathcg]:
            for i in range(0, 10):
                crossvalidationpath = os.path.join(base, str(i))
                os.mkdir(crossvalidationpath)
                forepath = os.path.join(crossvalidationpath, fore)
                os.mkdir(forepath)
                backpath = os.path.join(crossvalidationpath, back)
                os.mkdir(backpath)

        for _, _, y in os.walk('.'):
            all_files = y
            break

        for file in all_files:
            for i in range(0, 10):
                if 'foreground' in file and str(i) + '.' in file:
                    for base in [self.pathshuffle, self.pathcg]:
                        path1 = os.path.join(self.common_files, file)
                        path2 = os.path.join(base, str(i), fore, file)
                        shutil.copy(path1, path2)
                    os.remove(path1)
                elif 'background_cg' in file and str(i) + '.' in file:
                    path1 = os.path.join(self.common_files, file)
                    path2 = os.path.join(self.pathcg, str(i), back, file)
                    os.rename(path1, path2)
                elif 'background_biasaway' in file and str(i) + '.' in file:
                    path1 = os.path.join(self.common_files, file)
                    path2 = os.path.join(self.pathshuffle, str(i), back, file)
                    os.rename(path1, path2)

    def perform_TFFM_analysis(self, memetxt, backgroundpath, backgroundtype):

        script_name = os.path.join(self.script_base, 'TFFM_analysis.py')
        waitfile = os.path.join(self.common_files, 'waitTFFMfile.txt')
        meme = os.path.join(self.common_files, memetxt)

        for cross_validation in range(0, 10):
            session_name = self.species_name + self.tf_name + self.peaksize + str(cross_validation)

            folder_in_use = os.path.join(backgroundpath, str(cross_validation))

            training_fa_fore = os.path.join(folder_in_use, 'foreground',
                                            'training_' + 'foreground' + self.species_name + str(
                                                cross_validation) + '.fa')

            testing_fa_fore = os.path.join(folder_in_use, 'foreground',
                                           'testing_' + 'foreground' + self.species_name + str(
                                               cross_validation) + '.fa')

            testing_fa_back = os.path.join(folder_in_use, 'background',
                                           'testing_' + backgroundtype + self.species_name + str(
                                               cross_validation) + '.fa')

            with cd(folder_in_use):
                subprocess.call(['tmux', 'new-session', '-d', '-s', session_name])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', 'conda ', 'activate ',
                                 'tfbsshape', 'C-m'])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', 'python ', script_name, ' ',
                                 testing_fa_back, ' ',
                                 testing_fa_fore, ' ',
                                 training_fa_fore, ' ',
                                 meme,
                                 ' >>',
                                 'outcaptureTFFM.txt ', '2>&1 ', ' &&', ' echo ', session_name, ' >> ', waitfile, ' &&',
                                 ' tmux', ' kill-window', 'C-m'])

    def perform_PWM_analysis(self, memetxt, backgroundpath, backgroundtype):
        script_name = os.path.join(self.script_base, 'PWM_analysis.py')
        waitfile = os.path.join(self.common_files, 'waitPWMfile.txt')
        meme = os.path.join(self.common_files, memetxt)

        for cross_validation in range(0, 10):
            session_name = self.species_name + self.tf_name + self.peaksize + str(cross_validation)

            folder_in_use = os.path.join(backgroundpath, str(cross_validation))

            testing_fa_fore = os.path.join(folder_in_use, 'foreground',
                                           'testing_' + 'foreground' + self.species_name + str(
                                               cross_validation) + '.fa')

            testing_fa_back = os.path.join(folder_in_use, 'background',
                                           'testing_' + backgroundtype + self.species_name + str(
                                               cross_validation) + '.fa')

            with cd(folder_in_use):
                subprocess.call(['tmux', 'new-session', '-d', '-s', session_name])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', 'conda ', 'activate ',
                                 'tfbsshape', 'C-m'])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', 'python ', script_name, ' ',
                                 testing_fa_back, ' ',
                                 testing_fa_fore, ' ',
                                 meme,
                                 ' >>',
                                 'outcapturePWM.txt ', '2>&1 ', ' &&', ' echo ', session_name, ' >> ', waitfile, ' &&',
                                 ' tmux', ' kill-window', 'C-m'])

    def perform_shape_analysis(self, backgroundpath, backgroundtype):
        script_name = os.path.join(self.script_base, 'shape_analysis.sh')
        pssm = os.path.join(self.common_files, self.pfm)
        waitfile = os.path.join(self.common_files, 'waitshapeflie.txt')

        for cross_validation in range(0, 10):
            session_name = self.species_name + self.tf_name + self.peaksize + str(cross_validation)

            folder_in_use = os.path.join(backgroundpath, str(cross_validation))

            testing_fa_back = os.path.join(
                folder_in_use, 'background',
                'testing_' + backgroundtype + self.species_name + str(cross_validation) + '.fa')
            training_fa_back = os.path.join(
                folder_in_use, 'background',
                'training_' + backgroundtype + self.species_name + str(cross_validation) + '.fa')
            testing_bed_back = os.path.join(
                folder_in_use, 'background',
                'testing_' + backgroundtype + self.species_name + str(cross_validation) + '.bed')
            training_bed_back = os.path.join(
                folder_in_use, 'background',
                'training_' + backgroundtype + self.species_name + str(cross_validation) + '.bed')
            testing_fa_fore = os.path.join(
                folder_in_use, 'foreground',
                'testing_' + 'foreground' + self.species_name + str(cross_validation) + '.fa')
            training_fa_fore = os.path.join(
                folder_in_use, 'foreground',
                'training_' + 'foreground' + self.species_name + str(cross_validation) + '.fa')
            testing_bed_fore = os.path.join(
                folder_in_use, 'foreground',
                'testing_' + 'foreground' + self.species_name + str(cross_validation) + '.bed')
            training_bed_fore = os.path.join(
                folder_in_use, 'foreground',
                'training_' + 'foreground' + self.species_name + str(cross_validation) + '.bed')

            with cd(folder_in_use):
                tffm_f = 'tffm_first_order.xml'
                tffm_d = 'tffm_detailed.xml'
                subprocess.call(['tmux', 'new-session', '-d', '-s', session_name])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', 'conda ', 'activate ',
                                 'tfbsshape', 'C-m'])
                subprocess.call(['tmux', 'send-keys', '-t', session_name + ':0', script_name, ' ',
                                 testing_fa_back, ' ',
                                 training_fa_back, ' ',
                                 testing_bed_back, ' ',
                                 training_bed_back, ' ',
                                 testing_fa_fore, ' ',
                                 training_fa_fore, ' ',
                                 testing_bed_fore, ' ',
                                 training_bed_fore, ' ',
                                 pssm, ' ',
                                 tffm_f, ' ',
                                 tffm_d,
                                 ' >>',
                                 'outcaptureshape.txt ', '2>&1 ', ' &&', ' echo ', session_name, ' >> ', waitfile,
                                 ' &&',
                                 ' tmux', ' kill-window', 'C-m'])

    def waitabit(self, waitfile):
        waitfile = os.path.join(self.common_files, waitfile)
        tmux_working = True
        while tmux_working:
            try:
                with open(waitfile, 'r') as tmux_summary:
                    count = 0
                    for cv_number in tmux_summary:
                        count += 1
                    if count != 10:
                        time.sleep(5)
                        print('waiting')
                    else:
                        tmux_working = False
                        print('finished')
            except:
                time.sleep(5)
                print('waiting')

    def run(self):

        with cd("/home/bgp01/tfbsshape/data/clasification"):
            self.create_pfm_from_meme(self.memexml, self.pfm, self.tf_name, self.species_name)

        with cd("/home/bgp01/tfbsshape/data/"):
            self.create_folders()

        with cd(self.common_files):
            ###### generate foreground by split oiginal peaks filie
            ###### using python
            name = 'foreground' + self.species_name
            self.split_files_cv(self.true_peaks, 500, 10, name, bed=self.true_bed)

            ###### generate background by shuffling sequences in Ti or Fi while matching the same dinucleotide
            ###### composition using biasaway
            # self.generate_false_peaks_biasaway()
            # self.split_files(self.false_peaks_biasaway, 'background_biasaway',ç)

            ###### generate background  by randomly selecting the same number ofsequences as in Ti or Fi and
            ###### matching the same %GC composition distribution from a set of genomic background sequences
            ###### using python

            self.download_genome(self.url)
            self.calcularte_CG_content(self.true_peaks)
            self.generate_false_peaks_cg()
            name = 'background_cg' + self.species_name
            self.split_files_cv(self.randomfasta, 500, 10, name)

            ###### clasification of generated foreground and background files in folders for 10foldCV
            self.data_clasification()

            ###### clasification of generated foreground and background files in folders for 10foldCV
            self.perform_TFFM_analysis(self.memetxt, self.pathcg, 'background_cg')
            
            self.waitabit('waitTFFMfile.txt')
            
            self.perform_PWM_analysis(self.memetxt, self.pathcg, 'background_cg')
            
            self.waitabit('waitPWMfile.txt')
            
            self.perform_shape_analysis(self.pathcg, 'background_cg')

######### CAUTION, CHANGE SELF.GENOME IF FILE GENOME FILE IS NOT .GZ
shape_target = generate_training_files('Solanum_lycopersicum', 'MYC3', '51', 'SolycMYC3_51bp.correction.fa', 'SolycMYC3_51bp.bed',
                                       'meme.xml', 'meme.txt',
                                       'sol_shape',
                                       'ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/assembly/current_build/S_lycopersicum_chromosomes.4.00.fa.gz',
                                       '/home/bgp01/tfbsshape/scripts')

shape_target.run()
