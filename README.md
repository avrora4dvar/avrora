# AVRORA BALANCE

As used in Glider study. 

## AVRORA 4

1. Build the NETCDF3 library libnetcdf.a and place it in Avrora/Lib
2. Set MPI Fortran compiler and appropriate flags in the first block of the
   [makefile][Avrora/Compile/makefile]
3. cd into [Compile][Avrora/Compile]
4. Run './build.sh'

## ROMS
1. Copy [ak_or.h][Avrora/Include/ak_or.h] to 'ROMS/Include/' of the ROMS directory 
2. Compile the MPI, NETCDF 3 version of ROMS following instructions included in the ROMS distribution using ak_or.h as header file. 
3. Move the binary oceanM to [Bin][Avrora/Bin] as oceanMCR
4. Redirect the paths in the [ROMS input][Script/ocean.in] to the locations of the input files on the local system

** Scripts
1. Add scripts in [Script][Script] to [Bin][Avrora/Bin] and move to a local location for the Bin
2. Assign executive priviliges to Bin using `chmod`
3. Set bin directory (binDir) and experiment directory (expDir) in `batch_bal.sh`
4. Go through the othe scripts and set 
   -'tiles': the horizontal tiling for 'ad_avrora' and 'tl_avrora'
   -'np': number of threads used simultaneously per program. For 'tl_avrora' and 'ad_avrora' should be equal to product of 'tiles'
   -'nptot': total number of available cores for the calculations. Must be a multiple of 'np'.
5. Adjust the lines that start the different programs in the bin directory to the correct syntax for the system. 

## Experiment directory
1. Experiment directories names should be of the form <experiment name>_####
where #### is the start date of the window in number of days since 2005-01-01
2. Each experiment directory should contain a directory Iter0 containing the history files from the nonlinear in the previous window
3. Each experiment directory should contain a file 'case_def.txt' in which among others the paths to the boundary forcing files and the experiment number and start date should be set. See [Scripts][Scripts/case_def.txt] for an example.
4. Files with observations do not have to be stored in the experiment directory, their filename must be of the format obslist####.nc where #### is the start date of the associated window in days since 2005-01-01

## Run
batch_bal.sh ###1 ###2 > output.txt 
where ###1 is the number (in days since 2005-01-01) of the 1st window to be calcalculated and ###2 the last. ####2-####1 must be 0 mod 3. 

