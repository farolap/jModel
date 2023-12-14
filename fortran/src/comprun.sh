#!/bin/bash
set -e



export OMP_NUM_THREADS=15

FFLAGS='-O3 -fopenmp';
#FFLAGS=$FFLAGS' -fcheck=all -Wall';
grid='A';

test ! -d ../result${grid} && mkdir ../result${grid}

addit='';
if [[ $grid == "A" ]]; then
    swmop='swma_operators'
    timed='time_derivA'
    iniop='init_condA'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
    execut='shall_wat3dA'
elif [[ $grid == "Anorm" ]]; then
    swmop='swma_operators'
    timed='time_derivA'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj"
    execut='normmodesA'
elif [[ $grid == "B" ]]; then
    swmop='swmb_operators'
    iniop='init_condB'
    timed='time_derivB'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
    execut='shall_wat3dB'
elif [[ $grid == "Bnorm" ]]; then
    swmop='swmb_operators'
    iniop='init_condB'
    timed='time_derivB'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
    execut='normmodesB'
elif [[ $grid == "Ctri" ]]; then
    swmop='swmc_operators'
    iniop='init_condC'
    timed='time_derivC'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
    if [[ $grid == "Chex" ]]; then
        execut='shall_wat3dChex'
    elif [[ $grid == "Ctri" ]]; then
        execut='shall_wat3dCtri'
    fi
elif [[ $grid == "Ctrinorm" ]]; then
    swmop='swmc_operators'
    iniop='init_condC'
    timed='time_derivC'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj"
    execut='normmodesCtri'
elif [[ $grid == "Ainst" ]]; then
    swmop='swma_operators'
    iniop='init_condA'
    timed='time_derivA'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
	execut='instanaA'
elif [[ $grid == "Binst" ]]; then
    swmop='swmb_operators'
    iniop='init_condB'
    timed='time_derivB'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
	execut='instanaB'
elif [[ $grid == "Ctriinst" ]]; then
    swmop='swmc_operators'
    iniop='init_condC'
    timed='time_derivC'
    addit=" bkup/${swmop}.obj bkup/${timed}.obj bkup/${iniop}.obj bkup/init_cond.obj"
	execut='instanaCtri'
fi


# OBJ='bkup/constants.obj bkup/basic_funcs.obj bkup/datastruct.obj bkup/lmeshpack.obj bkup/smeshpack.obj ''bkup/'$swmop'.obj'


OBJ="bkup/constants.obj bkup/basic_funcs.obj bkup/datastruct.obj bkup/lmeshpack.obj bkup/smeshpack.obj${addit}"

gfortran $FFLAGS -march=native -c  constants.f90 -o constants.obj
gfortran $FFLAGS -march=native -c  datastruct.f90 -o datastruct.obj
gfortran $FFLAGS -march=native -c  basic_funcs.f90 -o basic_funcs.obj
gfortran $FFLAGS -march=native -c  lmeshpack.f90 -o lmeshpack.obj
gfortran $FFLAGS -march=native -c  smeshpack.f90 -o smeshpack.obj
if [[ $grid == "A" ]]; then
    gfortran $FFLAGS -march=native -c  swma_operators.f90 -o swma_operators.obj
    gfortran $FFLAGS -march=native -c  init_cond.f90 -o init_cond.obj
    gfortran $FFLAGS -march=native -c  init_condA.f90 -o init_condA.obj
    gfortran $FFLAGS -march=native -c  time_derivA.f90 -o time_derivA.obj
elif [[ $grid == "Anorm" ]]; then
    gfortran $FFLAGS -march=native -c  swma_operators.f90 -o swma_operators.obj
    gfortran $FFLAGS -march=native -c  time_derivA.f90 -o time_derivA.obj
elif [[ $grid == "B" ]]; then
    gfortran $FFLAGS -march=native -c  swmb_operators.f90 -o swmb_operators.obj
    gfortran $FFLAGS -march=native -c  init_cond.f90 -o init_cond.obj
    gfortran $FFLAGS -march=native -c  init_condB.f90 -o init_condB.obj
    gfortran $FFLAGS -march=native -c  time_derivB.f90 -o time_derivB.obj
elif [[ $grid == "Bnorm" ]]; then
    gfortran $FFLAGS -march=native -c  swmb_operators.f90 -o swmb_operators.obj
    gfortran $FFLAGS -march=native -c  time_derivB.f90 -o time_derivB.obj
elif [[ $grid == "Ctri" ]]; then
    gfortran $FFLAGS -march=native -c  swmc_operators.f90 -o swmc_operators.obj
    gfortran $FFLAGS -march=native -c  init_cond.f90 -o init_cond.obj
    gfortran $FFLAGS -march=native -c  init_condC.f90 -o init_condC.obj
    gfortran $FFLAGS -march=native -c  time_derivC.f90 -o time_derivC.obj
elif [[ $grid == "Ctrinorm" ]]; then
    gfortran $FFLAGS -march=native -c  swmc_operators.f90 -o swmc_operators.obj
    gfortran $FFLAGS -march=native -c  init_condC.f90 -o init_condC.obj
    gfortran $FFLAGS -march=native -c  time_derivC.f90 -o time_derivC.obj
elif [[ $grid == "Ainst" ]]; then
    gfortran $FFLAGS -march=native -c  swma_operators.f90 -o swma_operators.obj
    gfortran $FFLAGS -march=native -c  init_cond.f90 -o init_cond.obj
    gfortran $FFLAGS -march=native -c  init_condA.f90 -o init_condA.obj
    gfortran $FFLAGS -march=native -c  time_derivA.f90 -o time_derivA.obj
elif [[ $grid == "Binst" ]]; then
    gfortran $FFLAGS -march=native -c  swmb_operators.f90 -o swmb_operators.obj
    gfortran $FFLAGS -march=native -c  init_cond.f90 -o init_cond.obj
    gfortran $FFLAGS -march=native -c  init_condB.f90 -o init_condB.obj
    gfortran $FFLAGS -march=native -c  time_derivB.f90 -o time_derivB.obj
elif [[ $grid == "Ctriinst" ]]; then
    gfortran $FFLAGS -march=native -c  swmc_operators.f90 -o swmc_operators.obj
    gfortran $FFLAGS -march=native -c  init_condC.f90 -o init_condC.obj
    gfortran $FFLAGS -march=native -c  time_derivC.f90 -o time_derivC.obj
fi


mv *.obj ./bkup


if [[ $grid == "A" ]]; then
    gfortran $FFLAGS -march=native "${execut}.f90" $OBJ -o $execut
elif [[ $grid == "B" ]]; then
    gfortran $FFLAGS -march=native "${execut}.f90" $OBJ -o $execut
else
    gfortran $FFLAGS -march=native "${execut}.f90" $OBJ -o $execut
fi
#gfortran -O3 -fopenmp -march=native construct_grid.f90 $OBJ -o construct_grid
mv *.mod ./bkup


chmod +x ./$execut
./$execut
exit 1
