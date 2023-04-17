#!/bin/bash
set -e

#. ${HOME}/intel/oneapi/setvars.sh
ulimit -s unlimited

export OMP_NUM_THREADS=15

FFLAGS='-qopenmp -O3 -xhost';
grid='B';

swmop='swmb_operators'
timed='time_derivB'
addit=" bkup/${swmop}.obj bkup/${timed}.obj"
execut='shall_wat3dB'

OBJ="bkup/constants.obj bkup/basic_funcs.obj bkup/datastruct.obj bkup/lmeshpack.obj bkup/smeshpack.obj${addit}"

ifort $FFLAGS -c constants.f90 -o constants.obj 
ifort $FFLAGS -c  datastruct.f90 -o datastruct.obj
ifort $FFLAGS -c  basic_funcs.f90 -o basic_funcs.obj
ifort $FFLAGS -c  lmeshpack.f90 -o lmeshpack.obj
ifort $FFLAGS -c  smeshpack.f90 -o smeshpack.obj
ifort $FFLAGS -c  swmb_operators.f90 -o swmb_operators.obj
ifort $FFLAGS -c  time_derivB.f90 -o time_derivB.obj


mv *.obj ./bkup

ifort $FFLAGS "${execut}.f90" $OBJ -o $execut

mv *.mod ./bkup


chmod +x ./$execut
./$execut
exit 1
