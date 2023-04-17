#!/bin/bash
set -e


export OMP_NUM_THREADS=15

#FFLAGS='-O3 -fopenmp';
FFLAGS='-I/home/farolap/pkgs/usr/include/'
FIFLAGS='-fcheck=all -Wall -lnetcdff -lnetcdf'
LD_LIBRARY_PATH='/home/farolap/pkgs/usr/lib'

OBJ="bkup/constants.obj bkup/basic_funcs.obj bkup/datastruct.obj bkup/lmeshpack.obj bkup/smeshpack.obj"

gfortran $FIFLAGS -march=native -c  constants.f90 -o constants.obj
gfortran $FIFLAGS -march=native -c  datastruct.f90 -o datastruct.obj
gfortran $FIFLAGS -march=native -c  basic_funcs.f90 -o basic_funcs.obj
gfortran $FIFLAGS -march=native -c  lmeshpack.f90 -o lmeshpack.obj
gfortran $FIFLAGS -march=native -c  smeshpack.f90 -o smeshpack.obj
execut='save_gridnc'

mv *.obj ./bkup

gfortran  -o $execut -I/home/farolap/pkgs/usr/include/ "${execut}.f90" -L/home/farolap/pkgs/usr/lib/ $FIFLAGS -march=native $OBJ 


mv *.mod ./bkup


chmod +x ./$execut
./$execut
exit 1
