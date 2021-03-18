import sys

import TFFM.tffm_module as tffm_module
from TFFM.constants import TFFM_KIND

testing_fa_back = sys.argv[1]
testing_fa_fore = sys.argv[2]
memetxt = sys.argv[3]

pwm = tffm_module.tffm_from_meme(memetxt, TFFM_KIND.ZERO_ORDER)

print "1st-order best"
with open('foreground_Zo_class_PWM_bueno.tsv', "w") as outputfile:
    for hit in pwm.scan_sequences(testing_fa_fore, only_best=True):
        outputfile.write(str(hit) + "\n")

with open('background_Zo_class_PWM_bueno.tsv', "w") as outputfile:
    for hit in pwm.scan_sequences(testing_fa_back, only_best=True):
        outputfile.write(str(hit) + "\n")
