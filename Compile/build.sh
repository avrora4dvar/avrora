#!/bin/bash

binDir='../Bin'
cppFile=../Include/cppdefs.h
tmpFile=../Include/tmp.h
cp $cppFile $tmpFile

#Make Bin
if [ -d $binDir ]; then rmdir $binDir; fi
mkdir $binDir

#-----------------------------------------------------------------

echo 'Build TL_AVRORA'
echo '#define MPI' > $cppFile
echo '#undef ADJOINT' >> $cppFile
echo '#define TANGENT' >> $cppFile
echo '#undef NONLINEAR' >> $cppFile
cat $tmpFile >> $cppFile
make tl_avrora; make clean
make tl_interp
make tl_sample
make clean
mv tl_* $binDir

echo 'Build AD_AVRORA'
echo '#define MPI' > $cppFile
echo '#define ADJOINT' >> $cppFile
echo '#undef TANGENT' >> $cppFile
echo '#undef NONLINEAR' >> $cppFile
cat $tmpFile >> $cppFile
make ad_avrora; make clean
make ad_interp
make ad_sample
make clean
mv ad_* $binDir

echo 'Build TL_BALANCE'
echo '#undef MPI' > $cppFile
echo '#undef ADJOINT' >> $cppFile
echo '#define TANGENT' >> $cppFile
echo '#undef NONLINEAR' >> $cppFile
cat $tmpFile >> $cppFile
make tl_balance_3D; make clean
mv tl_balance_3D $binDir

echo 'Build AD_BALANCE'
echo '#undef MPI' > $cppFile
echo '#define ADJOINT' >> $cppFile
echo '#undef TANGENT' >> $cppFile
echo '#undef NONLINEAR' >> $cppFile
cat $tmpFile >> $cppFile
make ad_balance_3D; make clean
mv ad_balance_3D $binDir

echo 'Build other files'
echo '#define MPI' > $cppFile
echo '#undef ADJOINT' >> $cppFile
echo '#undef TANGENT' >> $cppFile
echo '#undef NONLINEAR' >> $cppFile
cat $tmpFile >> $cppFile
make cov_ini_temp_uni; make clean
mv cov_ini_temp_uni $binDir/cov_bal
make PCG; make clean; mv PCG $binDir/pcg
make create_background; mv create_background $binDir
make create_ini; mv create_ini $binDir
make create_pre; mv create_pre $binDir
make clean

cp $tmpFile $cppFile
rm $tmpFile