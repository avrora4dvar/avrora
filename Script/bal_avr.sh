#!/bin/bash   
                                                              
#SBATCH -J inner_loop                                                         
#SBATCH -o log_inner_%j                                                       
#SBATCH -e error_avr_%j
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=24                                       
#SBATCH -p compute                                                       
#SBATCH -t 06:00:00  

#back 8, ad_avr (6) 11,

#INNER_LOOP $1 $2 performs $iter_max inner loop iterations. The corrected
#initial condition file after iterations is written to ini.nc
#
#Input:
# $1: file name of the case_def file
# $2: number of the outer loop. To complete the outer loop run outer_loop.sh

source $1 
iter_last=( -1 -1 -1 -1 )

#tiles
npar=4
tiles=(2 3)
np=$(( ${tiles[0]}*${tiles[1]} ))
host=aruba2

iter_dir=$exp_dir/$(printf Iter%d $2)
prior_iter_dir=$exp_dir/$(printf Iter%d $(( $2-1 )) )
background_file=$prior_iter_dir/back.nc
ens_dir=$exp_dir/Ens

#machfile
mach_file=$exp_dir/mach.txt
echo $host > $mach_file

#--------------------------------------------------------------------------
#Sample forecast

if [ ${iter_last[0]} -lt 0 ]; then

 #Sample from nonlinear model for forecast                                    
 in_file=$prior_iter_dir/tl_sample_for.in
 out_file=$prior_iter_dir/tl_sample_for.out
 echo '#Observation list' > $in_file
 echo $obs_file >> $in_file
 echo '#History file' >> $in_file
 echo  $prior_iter_dir >> $in_file
 echo 'ocean_his' >> $in_file
 echo '#Grid file' >> $in_file
 echo  $roms_grid_file >> $in_file
 echo '#Input file' >> $in_file
 echo  $prior_iter_dir >> $in_file
 echo 'ocean_his' >> $in_file
 echo '#Output file' >> $in_file
 echo  $prior_iter_dir'/sample_for.nc' >> $in_file
 echo '#Output variable name' >> $in_file
 echo 'K' >> $in_file
 echo '#Reference time' >> $in_file
 echo $exp_time,$(( $exp_time+$window_duration )),$((-$exp_time)) >> $in_file

 if [ -f $out_file ]; then rm $out_file; fi
 echo start tl_sample at `date`
 ssh $host $bin_dir/tl_sample < $in_file > $out_file
 echo tl_sample done at `date`
 ln -s $prior_iter_dir/sample_for.nc $prior_iter_dir/sample_ana.nc

 #Create initial conditions next window                                        
 in_file=$prior_iter_dir'/ini.in'
 out_file=$prior_iter_dir'/ini.out'
 echo '#grid files' > $in_file
 echo $roms_grid_file >> $in_file
 echo $roms_grid_file >> $in_file
 echo $roms_grid_file >> $in_file
 echo '#netcdf files' >> $in_file
 echo $prior_iter_dir'/ocean_his_0003.nc' >> $in_file
 echo ''  >> $in_file
 echo $prior_iter_dir'/ini_for.nc' >> $in_file
 echo '#time steps' >> $in_file
 echo 24 1 >> $in_file
 echo '#weights' >> $in_file
 echo 1.0 0.0 >> $in_file
 echo '#output time' >> $in_file
 echo $(( $exp_time*86400 )) >> $in_file
 echo '#Default values' >> $in_file
 echo $T0,$S0,$Kt,$Km >> $in_file
 echo '#Roms runable' >> $in_file
 echo T >> $in_file

 if [ -f $out_file ]; then rm $out_file; fi
 echo start create_ini at `date`
 $bin_dir/create_ini < $in_file > $out_file
 echo create_ini  done at `date`
 ln -s $prior_iter_dir/ini_for.nc $prior_iter_dir/ini.nc

fi

#--------------------------------------------------------------------------
#Calculate corrections

work_dir=$iter_dir
def_file=$work_dir/case_def.txt

if [ ${iter_last[0]} -lt 0 ]; then

 #Create iteration directory
 if [ -d $iter_dir ]; then rm -r $iter_dir; fi
 mkdir $iter_dir

 #Create case_def file for this outer loop 
 cp $1 $def_file
 echo 'prior_ini_file='$prior_ini_file >> $def_file
 echo 'background_file='$background_file >> $def_file
 echo 'tiles=('${tiles[*]}')' >> $def_file
 echo 'ens_dir='$ens_dir >> $def_file

 #Copy samples from the background to the new outer loop directory
 ln -s $prior_iter_dir/sample_ana.nc $iter_dir/sample_for.nc
 cp $prior_iter_dir/sample_for.nc $iter_dir/r.nc
 ln -s $prior_iter_dir/ini.nc $iter_dir/ini_for.nc

 #Copy ensemble
 for iMem in `seq 1 999`; do
  member_dir=$exp_dir/Ens/$(printf Member_%03d $iMem )
  if [ -d $member_dir ]; then
   cp $member_dir/xobs.nc $member_dir/x.nc
   cp $member_dir/robs.nc $member_dir/r.nc
  fi
 done

 if [ ! -d $exp_dir/Pcg ]; then mkdir $exp_dir/Pcg; fi
 for iPar in `seq 1 $npar`; do
  pcg_dir=$exp_dir/Pcg/$(printf Member_%03d $iPar )
  if [ -d $pcg_dir ]; then rm -r $pcg_dir; fi
  mkdir $pcg_dir
 done
fi


#Iterations CGM
iter_max=14; iter=0 
while [ $iter -le $iter_max ]; do 

 echo STARTING INNER LOOP  $iter at `date`

 if [ $iter -gt ${iter_last[0]} ]; then

 #Set up administration for PCG
 in_file=$work_dir/pcg.in
 out_file=$work_dir/$(printf pcg%02d.out $iter)
 echo '#observation file' > $in_file
 echo $obs_file >> $in_file
 echo '#PCG file' >> $in_file
 echo $work_dir'/pcg.nc' >> $in_file
 echo 'File with B*residual' >> $in_file
 echo $work_dir >> $in_file
 echo $ens_dir >> $in_file
 echo $exp_dir/Pcg >> $in_file
 echo 'Read number thread' >> $in_file
 echo $npar >> $in_file
 echo 'Variable name in r_file' >> $in_file
 echo 'r.nc' >> $in_file
 echo 'K' >> $in_file
 echo 'Variable name x_file' >> $in_file
 echo 'x.nc' >> $in_file
 echo 'K' >> $in_file
 echo 'Factor normalized observations are marked as outliers' >> $in_file
 echo 10  >> $in_file
 echo 'Mode' >> $in_file
 echo 'B' >> $in_file

 if [ -f $out_file ]; then rm $out_file; fi
 echo 'running PCG'
 ssh $host $bin_dir/pcg < $in_file > $out_file
 fi

 #Exit loop 
 if [ $iter -eq $iter_max ]; then break; fi
# if [ -f $iter_dir/pcg_done.txt ]; then break; fi 

 #run ad_sample and ad_avrora
 if [ $iter -gt ${iter_last[1]} ]; then
  for iPar in `seq 1 $npar`; do
   pcg_dir=$exp_dir/Pcg/$(printf Member_%03d $iPar )
   pcg_no=$(( ($iPar-1)*$np ))
   cp $def_file $pcg_dir'/run_ad.in'
   echo 'input='$pcg_dir'/r.nc' >> $pcg_dir'/run_ad.in'
   echo 'output='$pcg_dir'/Ad.nc' >> $pcg_dir'/run_ad.in'
   $bin_dir/run_ad.sh $pcg_dir'/run_ad.in' $pcg_dir $host &
  done
  wait 
 fi

 #run covariance
 if [ $iter -gt ${iter_last[2]} ]; then
  for iPar in `seq 1 $npar`; do
   pcg_dir=$exp_dir/Pcg/$(printf Member_%03d $iPar )
   pcg_no=$(( ($iPar-1)*$np ))
   cp $def_file $pcg_dir'/run_cov.in'
   echo 'input='$pcg_dir'/Ad.nc' >> $pcg_dir'/run_cov.in'
   echo 'output='$pcg_dir'/Cb.nc' >> $pcg_dir'/run_cov.in'
   $bin_dir/run_cov_bal.sh $pcg_dir'/run_cov.in' $pcg_dir $host & 
  done
  wait
 fi

 #run tl_avrora and tl_sample
 if [ $iter -gt ${iter_last[3]} ]; then
  for iPar in `seq 1 $npar`; do
   pcg_dir=$exp_dir/Pcg/$(printf Member_%03d $iPar )
   pcg_no=$(( ($iPar-1)*$np ))
   cp $def_file $pcg_dir'/run_tl.in'
   echo 'input='$pcg_dir'/Cb.nc' >> $pcg_dir'/run_tl.in'
   echo 'output='$pcg_dir'/r.nc' >> $pcg_dir'/run_tl.in'
   $bin_dir/run_tl.sh $pcg_dir'/run_tl.in' $pcg_dir $host &
  done
  wait
 fi
 
 #update inner loop number
 iter=$(( $iter+1 ))
done 

#Remove CG
for iPar in `seq 1 $npar`; do
  pcg_dir=$exp_dir/Pcg/$(printf Member_%03d $iPar )
  #if [ -d $pcg_dir ]; then rm -r $pcg_dir; fi
done

#-------------------------------------------------------------------------
#Convert correction from dual to primal space

if [ $iter -gt ${iter_last[1]} ]; then
cp $def_file $iter_dir'/run_ad.in'
echo 'input='$iter_dir'/x.nc' >> $iter_dir'/run_ad.in'
echo 'output='$iter_dir'/Adx.nc' >> $iter_dir'/run_ad.in'
$bin_dir/run_ad.sh $iter_dir'/run_ad.in' $iter_dir $host 
fi 

if [ $iter -gt ${iter_last[2]} ]; then
cp $def_file $iter_dir'/run_cov.in'
echo 'input='$iter_dir'/Adx.nc' >> $iter_dir'/run_cov.in'
echo 'output='$iter_dir'/per_ana.nc' >> $iter_dir'/run_cov.in'
$bin_dir/run_cov_bal.sh $iter_dir'/run_cov.in' $iter_dir $host 
fi

if [ $iter -gt ${iter_last[3]} ]; then 
rm $iter_dir/Adx.nc

#Create initial conditions next window                                         
in_file=$iter_dir'/ini.in'
out_file=$iter_dir'/ini.out'
echo '#grid files' > $in_file
echo $roms_grid_file >> $in_file
echo $roms_grid_file >> $in_file
echo $roms_grid_file >> $in_file
echo '#netcdf files' >> $in_file
echo $iter_dir'/ini_for.nc' >> $in_file
echo $iter_dir'/per_ana.nc'  >> $in_file
echo $iter_dir'/ini.nc' >> $in_file
echo '#time steps' >> $in_file
echo 1 1 >> $in_file
echo '#weights' >> $in_file
echo 1.0 1.0 >> $in_file
echo '#output time' >> $in_file
echo $(( $exp_time*86400 )) >> $in_file
echo '#Default values' >> $in_file
echo $T0,$S0,$Kt,$Km >> $in_file
echo '#Roms runable' >> $in_file
echo T >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start create_ini at `date`
$bin_dir/create_ini < $in_file > $out_file
echo create_ini  done at `date`

fi 
#-------------------------------------------------------------------------
#Run ROMS

#Create ROMS input file                                                         
roms_in=$iter_dir'/ocean.in'
if [ -f $roms_in ]; then rm $roms_in; fi
DT=90.0
NFAST=30
DSTART=$exp_time
NTIMES=5760
NHIS=40; NDEFHIS=960
NAVG=960; NDEFAVG=2880; NTSAVG=0
tiles=(4 6)
np=$(( ${tiles[0]}*${tiles[1]} ))

#Copy ROMS file replacing the variables with their values                      
echo creating ROMS file $roms_in
while read -ra line; do
 for word in ${line[*]}; do
  if [[ $word == '$'* ]]; then
   word=$(eval echo $word)
  fi
  echo -n $word' ' >> $roms_in
 done
 echo '' >> $roms_in #newline                                                  
done < $roms_file

#Run roms                                                                     
run_dir=`pwd`
cd $iter_dir
echo start ROMS at `date`
roms_out=$iter_dir'/ocean.out'
if [ -f ocean.out ]; then rm ocean.out; fi                                     
echo 'ibrun:'$np
echo 'cd '$iter_dir > $iter_dir/run_roms.sh
echo 'export OMP_NUM_THREADS='$np >> $iter_dir/run_roms.sh 
echo $bin_dir/oceanOCR '<' $roms_in '>' $roms_out  >> $iter_dir/run_roms.sh
chmod u=rwx $iter_dir/run_roms.sh
ssh $host $iter_dir/run_roms.sh 
#mpirun -n $np -machinefile $mach_file $bin_dir/oceanMCR ocean.in > ocean.out
echo ROMS done at `date`
cd $run_dir

#------------------------------------------------------------------------
#Sample 

#Sample from nonlinear model for analysis                                     
in_file=$iter_dir/tl_sample_ana.in
out_file=$iter_dir/tl_sample_ana.out
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
echo  $iter_dir'/sample_ana.nc' >> $in_file
echo '#Output variable name' >> $in_file
echo 'K' >> $in_file
echo '#Reference time' >> $in_file
echo $exp_time,$(( $exp_time+$window_duration )),$((-$exp_time)) >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start tl_sample at `date`
ssh $host $bin_dir/tl_sample < $in_file > $out_file
echo tl_sample done at `date`

#-------------------------------------------------------------------------
#Create NL

if [ ! -d $exp_dir/NL ]; then mkdir $exp_dir/NL; fi 
mv $iter_dir/ocean* $exp_dir/NL

#-------------------------------------------------------------------------
#Create next window

next_exp_time=$(( $exp_time+$window_duration ))
next_exp_dir=${exp_dir//$exp_time/$next_exp_time}
next_obs_file=${obs_file//$exp_time/$next_exp_time}

#Create new window
mkdir $next_exp_dir
mkdir $next_exp_dir/Iter0

#Link ROMS
ln -s $exp_dir/NL/ocean_his_* $next_exp_dir/Iter0

#Create new case_def file                                                      
case_file=$next_exp_dir'/case_def.txt'
echo $case_file
if [ -f $case_file ]; then rm $case_file; fi
while read line; do
 if [[ $line  == 'exp_time='* ]]; then
  echo 'exp_time='$next_exp_time >> $case_file
 elif [[ $line == 'exp_dir='* ]]; then
  echo 'exp_dir='$next_exp_dir >> $case_file
 elif [[ $line == 'obs_file='* ]]; then
  echo 'obs_file='$next_obs_file >> $case_file
 else
  echo $line >> $case_file
 fi
done < $1


#--------------------------------------------------------------------------
#Close


echo window DONE at `date`
