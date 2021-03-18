import os
import numpy as np
import pandas
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import auc
import sys

global fold, leng
fold = f'/home/bgp01/tfbsshape/data/{str(sys.argv[1])}/'

class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_the_frame(filename, classlabel):
    with open(filename, 'r') as scoresfile:
        if '.txt' in filename:
            df = pandas.read_csv(scoresfile, sep='\s+', usecols=['proba'], header=0)
        elif '.tsv' in filename:
            df = pandas.read_csv(scoresfile, sep='\s+', usecols=[7], names=['proba'])

        #         print('original file shape =',df.shape)
        df.dropna().reset_index(drop=True)
        n_pred, _ = df.shape
        #         print('drop nan file shape =',df.shape)
        df['class'] = np.hstack(([classlabel] * n_pred))
        return df


def get_importantdata(filename, threashold):
    with open(filename, 'r') as scoresfile:
        if '.txt' in filename:
            df = pandas.read_csv(scoresfile, usecols=['peak_id', 'proba'], sep='\s+', header=0)
        elif '.tsv' in filename:
            df = pandas.read_csv(scoresfile, sep='\s+', usecols=['peak_id', 'proba'],
                                 names=['peak_id', 'start', 'end', 'strand', 'sequence', 'model', 'number', 'proba'])

        df.dropna().reset_index(drop=True)
        df = df.loc[df['proba'] >= threashold].reset_index(drop=True)
        return df


def get_bedfile(bed):
    df = pandas.read_csv(bed, sep="\t", names=['chr', 'start', 'end', 'peak_id'])

    return df


def get_folders():
    for a, x, y in os.walk('.'):
        return x, y


fig, ax = plt.subplots()

def plot_it(fore, back, name):
    subpath = '51/learning/cg/'
    common = '51/common_files/'
    with cd(fold + common):
        folders, files = get_folders()
        for fa in files:
            if '.bed' in fa and 'random' not in fa:
                bedname = fa
                print(bedname)
        beddf = get_bedfile(bedname)
        print('bedshape', beddf.shape)

    with cd(fold + subpath):
        folders, files = get_folders()

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 1000)
    for folder in folders:
        if not '000' in folder:
            with cd(fold + subpath + folder):
                dfroc = pandas.concat([
                    get_the_frame(fore, 1),
                    get_the_frame(back, 0)
                ], ignore_index=True)
            fpr, tpr, _ = roc_curve(dfroc['class'], dfroc['proba'])
            area = auc(fpr, tpr)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(area)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0

    for tpr, fpr in zip(mean_tpr, mean_fpr):
        if fpr <= 0.05:
            tprmax = tpr
            fprthreashold = fpr
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    title = str(fold.split('/')[-2])
    species = str(fold.split('/')[-3])
    condition = False
    for folder in folders:
        if not '000' in folder:
            with cd(fold + subpath + folder):
                if condition:
                    newdf = get_importantdata(fore, tprmax)
                    olddf = pandas.concat([olddf, newdf], ignore_index=True)

                else:
                    olddf = get_importantdata(fore, tprmax)
                    condition = True

    datafinal = pandas.merge(beddf, olddf, on="peak_id")
    print(datafinal, tprmax)
    with cd(fold + subpath):
        datafinal.to_csv(
            "{}_{}_{}_{}_{}.prueba.csv".format(
                species,
                title,
                name,
                'peaks',
                '05'
            ), sep='\t', header=False, index=False
        )
        with open('summary_result_all_data.bed', 'a+') as summary:
            summary.write(
                '{}\t{}\t{}\t{}\t{}\t{}\n'.format(species, title, name,  mean_auc, tprmax, fprthreashold))


plot_it('foreground_Zo_class_PWM_bueno.tsv', 'background_Zo_class_PWM_bueno.tsv', 'PWM')
plot_it('foreground_fo_class_TMMF.tsv', 'background_fo_class_TMMF.tsv', 'TFFM')
plot_it('DNAshapedPSSM_fg_predictions.txt', 'DNAshapedPSSM_bg_predictions.txt', 'PWM_+_shape')
plot_it('DNAshapedTFFM_fo_fg_predictions.txt', 'DNAshapedTFFM_fo_bg_predictions.txt', 'TFFM_+_shape')