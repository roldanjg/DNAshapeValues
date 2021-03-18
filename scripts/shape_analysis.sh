export PYTHONPATH=$PYTHONPATH:./
# Define where are located the DNA shape data from GBshape
araTha=/home/bgp01/tfbsshape/data/Solanum_lycopersicum/shape_vals/sol_shape
helt=$araTha/SL4.0.HelT.wig.bw;
mgw=$araTha/SL4.0.MGW.wig.bw;
prot=$araTha/SL4.0.ProT.wig.bw;
roll=$araTha/SL4.0.Roll.wig.bw;
helt2=$araTha/SL4.0.HelT.wig.bw;
mgw2=$araTha/SL4.0.MGW.wig.bw;
prot2=$araTha/SL4.0.ProT.wig.bw;
roll2=$araTha/SL4.0.Roll.wig.bw;
dnashaped=/home/bgp01/tfbsshape/programs/DNAshapedTFBS/DNAshapedTFBS.py;
testing_fa_back=$1
training_fa_back=$2
testing_bed_back=$3
training_bed_back=$4
testing_fa_fore=$5
training_fa_fore=$6
testing_bed_fore=$7
training_bed_fore=$8
pssm=$9
tffm_f=${10}
tffm_d=${11}

echo "Training a first order TFFM + DNA shape classifier.";
time python2.7 ${dnashaped} trainTFFM -T ${tffm_f} \
    -i ${training_fa_fore} \
    -I ${training_bed_fore} \
    -b ${training_fa_back} \
    -B ${training_bed_back} \
    -o DNAshapedTFFM_fo_classifier -t first_order \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
    
# echo "Training a detailed TFFM + DNA shape classifier.";
# time python2.7 ${dnashaped} trainTFFM -T ${tffm_d} \
#     -i ${training_fa_fore} -I ${training_bed_fore} \
#     -b ${training_fa_back} -B ${training_bed_back} \
#     -o DNAshapedTFFM_d_classifier -t detailed \
#     -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
    
echo "Training a PSSM + DNA shape classifier.";
time python2.7 ${dnashaped} trainPSSM -f ${pssm} \
    -i ${training_fa_fore} \
    -I ${training_bed_fore} \
    -b ${training_fa_back} \
    -B ${training_bed_back} \
    -o DNAshapedPSSM_classifier \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;

# echo "Training a 4-bits + DNA shape classifier.";
# time python2.7 ${dnashaped} train4bits -f ${pssm} \
#     -i ${training_fa_fore} \
#     -I ${training_bed_fore} \
#     -b ${training_fa_back} \
#     -B ${training_bed_back} \
#     -o DNAshaped4bits_classifier \
#     -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;

echo "Applying the trained first order TFFM + DNA shape classifier on foreground sequences.";
time python2.7 ${dnashaped} applyTFFM -T ${tffm_f} \
    -i ${testing_fa_fore} \
    -I ${testing_bed_fore} \
    -c DNAshapedTFFM_fo_classifier.pkl -o DNAshapedTFFM_fo_fg_predictions.txt \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

#echo "Applying the trained detailed TFFM + DNA shape classifier on foreground sequences.";
#time python2.7 ${dnashaped} applyTFFM -T ${tffm_d} \
#    -i ${testing_fa_fore} -I ${testing_bed_fore} \
#    -c DNAshapedTFFM_d_classifier.pkl -o DNAshapedTFFM_d_fg_predictions.txt \
#    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

echo "Applying the trained PSSM + DNA shape classifier on foreground sequences.";
time python2.7 ${dnashaped} applyPSSM -f ${pssm} \
    -i ${testing_fa_fore} \
    -I ${testing_bed_fore} \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_fg_predictions.txt \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

# echo "Applying the trained 4-bits + DNA shape classifier on foreground sequences.";
# time python2.7 ${dnashaped} apply4bits -f ${pssm} \
#     -i ${testing_fa_fore} \
#     -I ${testing_bed_fore} \
#     -c DNAshaped4bits_classifier.pkl -o DNAshaped4bits_fg_predictions.txt \
#     -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

echo "Applying the trained first order TFFM + DNA shape classifier on background sequences.";
time python2.7 ${dnashaped} applyTFFM -T ${tffm_f} \
    -i ${testing_fa_back} \
    -I ${testing_bed_back} \
    -c DNAshapedTFFM_fo_classifier.pkl -o DNAshapedTFFM_fo_bg_predictions.txt \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

#echo "Applying the trained detailed TFFM + DNA shape classifier on background sequences.";
#time python2.7 ${dnashaped} applyTFFM -T ${tffm_d} \
#    -i ${testing_fa_back} -I ${testing_bed_back} \
#    -c DNAshapedTFFM_d_classifier.pkl -o DNAshapedTFFM_d_bg_predictions.txt \
#    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;
    
echo "Applying the trained PSSM + DNA shape classifier on background sequences.";
time python2.7 ${dnashaped} applyPSSM -f ${pssm} \
    -i ${testing_fa_back} \
    -I ${testing_bed_back} \
    -c DNAshapedPSSM_classifier.pkl -o DNAshapedPSSM_bg_predictions.txt \
    -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;

# echo "Applying the trained 4-bits + DNA shape classifier on background sequences.";
# time python2.7 ${dnashaped} apply4bits -f ${pssm} \
#     -i ${testing_fa_back} \
#     -I ${testing_bed_back} \
#     -c DNAshaped4bits_classifier.pkl -o DNAshaped4bits_bg_predictions.txt \
#     -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0;
