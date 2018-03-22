#!/bin/bash

#SBATCH -J member                                                          
#SBATCH -o log_member_%j                                                    
#SBATCH --ntasks-per-node=6                                                   
#SBATCH --nodes=1                                                            
#SBATCH -p shared                                                             
#SBATCH -t 03:45:00                                                           
#SBATCH --export=ALL           
#SBATCH -A osu116 

#Convert 4DVAR correction from dual to primal space and move  member forward in time
#$1: case_def.txt file from last outer loop iteration
#$2: ensemble number
#$3: host

#Import definitions
source $1
tiles=(2 3)
np=$(( ${tiles[0]}*${tiles[1]} ))

echo 'STARTING MEMBER '$2

#Find latest iteration
for iter in `seq 9 -1 0`; do
 iter_dir=$exp_dir/Iter$iter
 if [ -d $iter_dir ]; then break; fi
done
iter_dir=$exp_dir/Iter$iter
echo 'Last iteration in '$iter_dir

ens_dir=$exp_dir/Ens
if [ ! -d $ens_dir ]; then mkdir $ens_dir; fi
work_dir=$(printf '%s/Member_%-.3d/' $ens_dir $2)
if [ ! -d $work_dir ]; then mkdir $work_dir; fi
echo 'Work directory:'$work_dir

#Write case_def file
def_file=$work_dir/case_def.txt
cp $1 $def_file
echo 'tiles=('${tiles[0]}' '${tiles[1]}')' >> $def_file
echo 'background_file='$iter_dir/back.nc >> $def_file
echo 'member_file=per_for.nc' >> $def_file
echo 'ens_dir='$ens_dir >> $def_file

#-----------------------------------------------------------------------------

#Apply adjoint 
in_file=$work_dir/run_ad.in
cp $def_file $in_file
echo 'input='$work_dir'/x.nc' >> $in_file
echo 'output='$work_dir'/Ad.nc' >> $in_file
$bin_dir/run_ad.sh $in_file $work_dir $3

#Apply covariance 
in_file=$work_dir/run_cov.in
cp $def_file $in_file
echo 'input='$work_dir'/Ad.nc' >> $in_file
echo 'output='$work_dir'/Cb.nc' >> $in_file
$bin_dir'/run_cov_bal.sh' $in_file $work_dir $3

#Apply tangent linear
in_file=$work_dir'/run_tl.in'
cp $def_file $in_file
echo 'input='$work_dir'/Cb.nc' >> $in_file
echo 'output='$work_dir'/robs.nc' >> $in_file
$bin_dir/run_tl.sh $in_file $work_dir $3 

#Clean
rm -f $work_dir/Ad.nc
rm -f $work_dir/Cb.nc

echo 'MEMBER '$2' DONE'
