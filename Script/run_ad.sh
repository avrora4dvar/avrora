#Runs ad_sample, ad_avrora and interpolation
#$1: case_def file with 
#$2: work directory
#$3: machine on which ad_avrora is executed (e.g. aruba2)

source $1 #case_def file
def_file=$1
work_dir=$2
host=$3
np=$(( ${tiles[0]} * ${tiles[1]} ))

#mach file
mach_file=$work_dir/mach.txt
echo $host > $mach_file

#Convert observation space -> model space (ad_sample)                          
in_file=$work_dir'/ad_sample.in'
out_file=$work_dir'/ad_sample.out'
echo '#Observation list' > $in_file
echo $obs_file >> $in_file
echo '#Background file' >> $in_file
echo $background_file >> $in_file
echo '#Grid file' >> $in_file
echo $avr_grid_file >> $in_file
echo '#Input file (coefficient of impulses at observation locations)' >> $in_file
echo $input >> $in_file
echo '#Variable name in output' >> $in_file
echo 'K' >> $in_file
echo '#Output file' >> $in_file
echo $work_dir'/as.nc' >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start ad_sample at `date`
ssh $host $bin_dir/ad_sample < $in_file > $out_file
echo ad_sample done at `date`

#ad_avrora
in_file=$work_dir'/ad_avr.in'
out_file=$work_dir'/ad_avr.out'

if [ -f $work_dir'/ai.nc' ]; then rm $work_dir'/ai.nc'; fi
if [ -f $work_dir'/af.nc' ]; then rm $work_dir'/af.nc'; fi
cp $prior_frc_file $work_dir'/af.nc'

echo '# Number of b/c time steps, NTIMES:' >$in_file
echo $NTIMES >> $in_file
echo 'B/c time step DT' >> $in_file
echo $DT >> $in_file
echo '# NFAST' >> $in_file
echo $NFAST >> $in_file
echo '# Output to history file each NHIS times' >> $in_file
echo $NHIS >> $in_file
echo '# Background eddy viscosity' >> $in_file
echo $Km >> $in_file
echo '# Background eddy diffusivity' >> $in_file
echo $Kt $Kt >> $in_file
echo '# Bottom friction coefficient' >> $in_file
echo $RDRG >> $in_file
echo '# Horizontal viscosity' >> $in_file
echo $VISC2 >> $in_file
echo '# Horizontal diffusion' >> $in_file
echo $DIFF2 >> $in_file
echo '# Params of linear EOS' >> $in_file
echo $R0 $TCOEF $SCOEF $T0 $S0 >> $in_file
echo '# Grid name' >> $in_file
echo $avr_grid_file >> $in_file
echo '# Initial condition file' >> $in_file
echo $work_dir'/ai.nc' >> $in_file
echo '#  Forcing file name' >> $in_file
echo $work_dir'/af.nc' >> $in_file
echo '# FWD (background) state:' >> $in_file
echo $background_file >> $in_file
echo '# History file' >> $in_file
echo $work_dir'/as.nc' >> $in_file
echo '# tiles' >> $in_file
echo ${tiles[*]} >> $in_file
#echo '1 1' >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start ad_avrora at `date`
mpirun -np $np -machinefile $mach_file $bin_dir'/ad_avrora' < $in_file > $out_file
echo ad_avrora done at `date`

#Interpolate 
in_file=$work_dir'/ad_interp.in'
out_file=$work_dir'/ad_interp.out'                                            
echo '#grid files' > $in_file
echo $avr_grid_file >> $in_file
echo $roms_grid_file >> $in_file
echo '#netcdf files' >> $in_file
echo $work_dir'/ai.nc' >> $in_file
echo $output >> $in_file
echo '#time steps' >> $in_file
echo 1 1 >> $in_file
echo '#output time' >> $in_file
echo 0.0 >> $in_file
echo '#Background ' >> $in_file
echo $background_file >> $in_file

if [ -f $out_file ]; then rm $out_file; fi
echo start ad_interp at `date`
ssh $host $bin_dir/ad_interp < $in_file > $out_file
echo ad_interp done at `date`

rm -f $work_dir/as.nc
rm -f $work_dir/ai.nc
rm -f $work_dir/af.nc






