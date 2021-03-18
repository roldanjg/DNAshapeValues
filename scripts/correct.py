import os
import subprocess
import time
class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def perform_PWM_analysis( memetxt, backgroundtype,species_name, tf):
    script_base = '/home/bgp01/tfbsshape/scripts'
    backgroundpath ='/home/bgp01/tfbsshape/data/'+species_name +'/' +tf + '/51/learning/cg'
    common_files = '/home/bgp01/tfbsshape/data/'+species_name +'/' +tf + '/51/common_files'
    script_name = os.path.join(script_base, 'PWM_analysis.py')
    waitfile = os.path.join(common_files, 'waitPWMBUENAfile.txt')
    meme = os.path.join(common_files, memetxt)

    for cross_validation in range(0, 10):
        session_name = species_name + tf + '51' + str(cross_validation)

        folder_in_use = os.path.join(backgroundpath, str(cross_validation))

        testing_fa_fore = os.path.join(folder_in_use, 'foreground',
                                       'testing_' + 'foreground' + species_name + str(
                                           cross_validation) + '.fa')

        testing_fa_back = os.path.join(folder_in_use, 'background',
                                       'testing_' + backgroundtype + species_name + str(
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


def waitabit(waitfile, species_name, tf):
    common_files = '/home/bgp01/tfbsshape/data/' + species_name + '/' + tf + '/51/common_files'
    waitfile = os.path.join(common_files, waitfile)
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






listtrans = ['bH018',  'bH031',  'bH066',  'BIM2',  'MYC2',  'MYC2CHLP',  'MYC3',  'MYC4',  'MYCH7',  'PIF4',  'PIF5']
for trans in listtrans:

    perform_PWM_analysis('meme.txt', 'background_cg','Arabidopsis_thaliana', trans)
    waitabit('waitPWMBUENAfile.txt','Arabidopsis_thaliana',trans)