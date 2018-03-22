#!/bin/bash
#RUN_COV_BAL runs different components of balance operator
#$1: case_def file
#$2: work directory which holds input/output files and tmp files
#$3: number host cpu

source $1
def_file=$1
work_dir=$2
host=$3

#mach file
mach_file=$work_dir/mach.txt
echo $host > $mach_file

#Apply ad_balance_3D                                                           
in_file=$work_dir/cov_ad.in
out_file=$work_dir/cov_ad.out
echo 'grid file' > $in_file
echo $roms_grid_file >> $in_file
echo 'Output file with temp,psi' >> $in_file
echo $work_dir'/tmp_Abal.nc' >> $in_file
echo 'Input file with temp,salt,u,v,zeta' >> $in_file
echo $input >> $in_file
echo 'Background file' >> $in_file
echo $background_file >> $in_file
echo 'R0,Tcoef,Scoef' >> $in_file
echo $R0,$TCOEF,$SCOEF >> $in_file

echo start cov_bal at `date` 
if [ -f $out_file ]; then rm $out_file; fi
ssh $host $bin_dir/ad_balance_3D < $in_file > $out_file
echo ad_balance_3D done at `date`

#Apply cov_bal                                                                
in_file=$work_dir'/cov_bal.in'
out_file=$work_dir'/cov_bal.out'
cp $work_dir'/tmp_Abal.nc' $work_dir'/tmp_Cb.nc'                                          
echo 'grid file' >  $in_file
echo $roms_grid_file >> $in_file
echo 'Output file with temp,psi' >> $in_file
echo $work_dir'/tmp_Cb.nc' >> $in_file
echo 'dlx' >> $in_file
echo $dlx >> $in_file
echo 'dly' >> $in_file
echo $dly >> $in_file
echo 'dlz' >> $in_file
echo $dlz >> $in_file
echo 'standard deviation sea-surface temperature errors' >> $in_file
echo $sig_file >> $in_file

echo start cov_bal at `date`
if [ -f $out_file ]; then rm $out_file; fi
ssh $host $bin_dir/cov_bal < $in_file > $out_file
echo cov_bal done at `date`

#Apply tl_balance_3D                                                           
in_file=$work_dir'/cov_tl.in'
out_file=$work_dir'/cov_tl.out'
echo 'grid file' > $in_file
echo $roms_grid_file >> $in_file
echo 'Input file with temp,salt,psi' >> $in_file
echo $work_dir'/tmp_Cb.nc' >> $in_file
echo 'Output file with temp,salt,u,v,zeta' >> $in_file
echo $output >> $in_file
echo 'Background file (not necessary)' >> $in_file
echo '' >> $in_file
echo 'R0,Tcoef,Scoef' >> $in_file
echo $R0,$TCOEF,$SCOEF >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start tl_balance_3D at `date`
ssh $host $bin_dir/tl_balance_3D < $in_file > $out_file
echo tl_balance_3D done at `date`

rm -f $work_dir/tmp*
