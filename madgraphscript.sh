#!/usr/bin/env bash
#
#THIS SCRIPT CREATES A DIRECTORY AND THEN RUNS MADGRAPH A NUMBER OF TIMES DEPENDING ON HOW MANY
#DIFFERENT EPSILON VALUES AND MZ VALUES YOU WANT TO TEST, EACH RUN OF MADGRAPH IS STORED IN
#THE INITIAL DIRECTORY CREATED
#

OUTPUT_DIR="scan_results" #CAN CHANGE THIS NAME TO WHATEVER YOWU WANT
mkdir -p "$OUTPUT_DIR"
mzMASSES=(1.0 2.0) #DIFFERENT MZ MASSES TO TEST
eMASSES=(1e-3) #DIFFERENT EPSILONS TO TEST
for mzMASS in "${mzMASSES[@]}"
do
	for eMASS in "${eMASSES[@]}"
	do
		RUNAME="mz_${mzMASS}_e_${eMASS}"
#
#THIS IS THE MADGRAPH PROCESS RUNNING
#
./mg5_aMC <<EOF
import model HAHM_variablesw_v3_UFO
generate p e- > p e- zp, zp > e+ e-
output ${RUNAME}
launch
1
0
#
#THIS SETS THE RUN CARD
#
set nevents 100
set ebeam1 120.0
set ebeam2 18.0
set lpp1 1
set lpp2 0
set time_of_flight 0.0
#
#THIS SETS THE PARAM CARD
#
set param_card mzdinput ${mzMASS} #mzdinput
set param_card DECAY 1023 Auto
set param_card HIDDEN 1 ${mzMASS} #mZDinput
set param_card HIDDEN 3 ${eMASS} #epsilon
0
exit
EOF
mv "$RUNAME" "$OUTPUT_DIR/"
done
done
