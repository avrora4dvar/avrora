#!/bin/bash

bin_dir=/home/aruba/vol2/ipasmans/Exp/Bin_Exp35

for iMem in `seq $1 3 $2`; do
    exp_dir=/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_$iMem
    def_file=$exp_dir/case_def.txt
    cd $exp_dir
    echo 'Ensemble '$iMem
    $bin_dir/batch_bal_member.sh $def_file 1 8
    echo 'PCG '$iMem
    $bin_dir/bal_avr.sh $def_file 1   
done