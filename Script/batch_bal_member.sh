#!/bin/bash

#SBATCH -J memberBatch                                                    
#SBATCH -o log_memberBatch_%j                                                 
#SBATCH --ntasks-per-node=12                                                  
#SBATCH --nodes=1                                                              
#SBATCH -p shared                                                            
#SBATCH -t 03:40:00                                                         
#SBATCH --export=ALL                                                        
#SBATCH -A osu116

#$1: case_def file
#$2: first representer to be calculated
#$3: last representer to be calculated

#Read input
source $1
host=aruba2
nptot=24; np=6
flag_copy=0

echo 'start batch_avr_member at '`date`

#Machfile
mach_file=$exp_dir/mach.txt
echo $host > $mach_file

#Find last iteration
for nIter in `seq 9 -1 0`; do
 iter_dir=$exp_dir/$(printf Iter%d $nIter)
 if [ -d $iter_dir ]; then break; fi
done
echo 'Iter0 dir: '$iter_dir

#------------------------------------------------------------------------
# Create background

background_file=$iter_dir/back.nc

#Calculate background                                                          
in_file=$iter_dir'/back.in'
out_file=$iter_dir'/back.out'
echo '#Grid file' > $in_file
echo $roms_grid_file >> $in_file
echo $avr_grid_file >> $in_file
echo '#Directory with history files' >> $in_file
echo $iter_dir >> $in_file
echo '#Number of first and last history file to be read' >> $in_file
echo 3 6 >> $in_file
echo '#Output background file' >> $in_file
echo $background_file >> $in_file
echo '#First time to be written to output (days since dateref)' >> $in_file
echo $exp_time $(( $exp_time + $window_duration )) 3600 >> $in_file
echo '#Reference time history files (yyyy mm dd HH MM SS)' >> $in_file
echo 2005 1 1 0 0 0 >> $in_file
echo '#Shift in reference time in output (days)' >> $in_file
echo $exp_time >> $in_file
echo '#Background diffusivity' >> $in_file
echo $T0,$S0,$Kt, $Km >> $in_file
echo '#Flag to activate filtering (T/F)' >> $in_file
echo F >> $in_file
echo '#Tiles' >> $in_file
echo 1,1 >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start create_background at `date`
$bin_dir/create_background < $in_file > $out_file
echo create_background done at `date`

#----------------------------------------------------------------------------
#Sample

#Sample from nonlinear model for analysis                                     
in_file=$iter_dir/tl_sample_for.in
out_file=$iter_dir/tl_sample_for.out
 echo '#Observation list' > $in_file
 echo $obs_file >> $in_file
 echo '#History file' >> $in_file
 echo  $iter_dir >> $in_file
 echo 'ocean_his' >> $in_file
 echo '#Grid file' >> $in_file
 echo  $roms_grid_file >> $in_file
 echo '#Input file' >> $in_file
 echo  $iter_dir >> $in_file
 echo 'ocean_his' >> $in_file
 echo '#Output file' >> $in_file
 echo  $iter_dir'/sample_for.nc' >> $in_file
 echo '#Output variable name' >> $in_file
 echo 'K' >> $in_file
 echo '#Reference time' >> $in_file
 echo $exp_time,$(( $exp_time+$window_duration )),$((-$exp_time)) >> $in_file

 if [ -f $out_file ]; then rm $out_file; fi
 echo start tl_sample at `date`
 $bin_dir/tl_sample < $in_file > $out_file
 echo tl_sample done at `date`

 ln -s $iter_dir/sample_for.nc $iter_dir/sample_ana.nc
 cp $iter_dir/sample_for.nc $iter_dir/r.nc

#---------------------------------------------------------------------------
#Preconditioner

ens_dir=$exp_dir/Ens
if [ ! -d $ens_dir ]; then
 mkdir $ens_dir

 for iMem in `seq 1 $nMembers`; do
  member_dir=$ens_dir/$(printf Member_%03d $iMem)
  mkdir $member_dir
 done
fi

#Create normal distributed vectors                                            
in_file=$ens_dir'/pre.in'
out_file=$ens_dir'/pre.out'
echo '#observation file' > $in_file
echo $obs_file >> $in_file
echo '#Input directories' >> $in_file
echo $iter_dir >> $in_file
echo $ens_dir >> $in_file
echo '#Read number thread' >> $in_file
echo $nMembers >> $in_file
echo '#Variable name in r_file' >> $in_file
echo 'r.nc' >> $in_file
echo 'K' >> $in_file
echo '#Variable name x_file' >> $in_file
echo 'xobs.nc' >> $in_file
echo 'K' >> $in_file
echo 'F#actor normalized observations are marked as outliers' >> $in_file
echo  10 >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start create_pre at `date`
ssh $host $bin_dir/create_pre < $in_file > $out_file
echo create_pre done at `date`

for iMem in `seq 1 $nMembers`; do
 member_dir=$exp_dir/Ens/$(printf Member_%03d $iMem)
 cp $member_dir/xobs.nc $member_dir/x.nc
done

#----------------------------------------------------------------------------
#Copy to local storage

if [ $flag_copy -ne 0 ]; then

 #Create exp_dir on local disk
 store_exp_dir=$exp_dir
 local_exp_dir=${exp_dir//$store_dir/$local_dir}
 echo 'creating '$local_exp_dir
 mkdir $local_exp_dir
 mkdir $local_exp_dir/Ens

 #Replace case_def
 echo creating case_def file
 def_file=$local_exp_dir/case_def.txt
 if [ -f $def_file ]; then rm $def_file; fi 
 while read -r line; do
  if [[ $line == 'exp_dir='* ]]; then
   echo ${line//$store_dir/$local_dir} >> $def_file
  elif [[ $line == 'ens_dir='* ]]; then 
   echo ${line//$store_dir/$local_dir} >> $def_file 
  elif [[ $line == 'background_file='* ]]; then 
   echo ${line//$store_dir/$local_dir} >> $def_file
  elif [[ $line == 'prior_ini_file='* ]]; then 
   echo ${line//$store_dir/$local_dir} >> $def_file
  else
   echo $line >> $def_file
  fi                                              
 done < $1

 cp $local_exp_dir/case_def.txt $store_exp_dir/case_local.txt

 #Copy members
 for iMem in `seq 1 $nMembers`; do
  prior_member_dir=$store_exp_dir/Ens/$(printf Member_%03d $iMem)
  member_dir=${prior_member_dir//$store_dir/$local_dir}
  echo 'copying '$prior_member_dir
  cp -r $prior_member_dir $member_dir
 done

 #Copy iter dir
 prior_member_dir=$iter_dir
 member_dir=${prior_member_dir//$store_dir/$local_dir}
 echo 'copying '$prior_member_dir
 cp -r $prior_member_dir $member_dir

else

 def_file=$1
 
fi #flag_copy


#--------------------------------------------------------------------------
#Start members

npactive=0
for iMem in `seq $2 $3`; do

 $bin_dir/bal_avr_member.sh $def_file $iMem $host &
 npactive=$(( $npactive+$np ))
 sleep 2

 if [ $npactive -eq $nptot ]; then
  wait
  npactive=0
 fi

done
wait

#------------------------------------------------------------------------
#Copy back

if [ $flag_copy -ne 0 ]; then

 for iMem in `seq $2 $3`; do
 
  prior_member_dir=$local_exp_dir/Ens/$(printf Member_%03d $iMem)
  member_dir=$store_exp_dir/Ens/$(printf Member_%03d $iMem)
  echo 'Copying back '$member_dir
  mv $prior_member_dir/* $member_dir
 done

fi #flag_copy

#------------------------------------------------------------------------

echo 'batch_bal_member DONE at '`date`
