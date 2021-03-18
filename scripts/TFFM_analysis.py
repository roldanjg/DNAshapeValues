import sys
import TFFM.tffm_module as tffm_module
from TFFM.constants import TFFM_KIND

testing_fa_back = sys.argv[1]
testing_fa_fore = sys.argv[2]
training_fa_fore = sys.argv[3]
print training_fa_fore
memetxt = sys.argv[4]

print "1st-order best"
tffm_first_order = tffm_module.tffm_from_meme(memetxt, TFFM_KIND.FIRST_ORDER)
tffm_first_order.write("tffm_first_order_initial.xml")
print "1st-order best"
tffm_first_order.train(training_fa_fore)
tffm_first_order.write("tffm_first_order.xml")

tffm_detailed = tffm_module.tffm_from_meme(memetxt, TFFM_KIND.DETAILED)
tffm_detailed.write("tffm_detailed_initial.xml")
print "1st-order best"
tffm_detailed.train(training_fa_fore)
print "1st-order best"
tffm_detailed.write("tffm_detailed.xml")


tffm_first_order = tffm_module.tffm_from_xml("tffm_first_order.xml",
        TFFM_KIND.FIRST_ORDER)

print "1st-order best"
with open('foreground_fo_class_TMMF.tsv', "w") as outputfile:
    for hit in tffm_first_order.scan_sequences(testing_fa_fore, only_best=True):
        outputfile.write(str(hit) + "\n")

with open('background_fo_class_TMMF.tsv', "w") as outputfile:
    for hit in tffm_first_order.scan_sequences(testing_fa_back, only_best=True):
        outputfile.write(str(hit) + "\n")

tffm_detailed = tffm_module.tffm_from_xml("tffm_detailed.xml",
        TFFM_KIND.DETAILED)

print "detailed best"
with open('foreground_fo_class_TMMF.tsv', "w") as outputfile:
    for hit in tffm_detailed.scan_sequences(testing_fa_fore, only_best=True):
        outputfile.write(str(hit) + "\n")

with open('background_fo_class_TMMF.tsv', "w") as outputfile:
    for hit in tffm_detailed.scan_sequences(testing_fa_back, only_best=True):
        outputfile.write(str(hit) + "\n")