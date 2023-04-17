#----------------------------------
# Makefile
#----------------------------------

OS := $(shell uname -s)
F90:= $(shell which gfortran)
FFLAG := -O3 -fopenmp -march=native
model:= construct_grid.f90

IMOD=-Ibin
#Objects
OBJ =   bin/constants.obj \
		bin/datastruct.obj \
		bin/basic_funcs.obj \
		bin/lmeshpack.obj \
		bin/smeshpack.obj \



all: header config bin/shw  ending

run: all
	./shw

header:
	@echo --------------------------------------------
	@echo Compiling and building the software   
	@echo --------------------------------------------

#Configure Enviroment (directories)
config:
	chmod +x src/*.sh
	. src/dirs.sh

#Constants
bin/constants.obj: src/constants.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv constants.mod bin/.

bin/datastruct.obj: src/datastruct.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv datastruct.mod bin/.

bin/basic_funcs.obj: src/basic_funcs.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv basic_funcs.mod bin/.

bin/lmeshpack.obj: src/lmeshpack.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv lmeshpack.mod bin/.

bin/smeshpack.obj: src/smeshpack.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv smeshpack.mod bin/.

# bin/swma_operators.obj: src/swma_operators.f90
# 	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
# 	mv swma_operators.mod bin/.


#Main executable
bin/shw: src/$(model) $(OBJ)
	$(F90) $(FFLAG)  src/$(model) $(OBJ) -o $@ $(IMOD)

ending: 
	chmod +x src/link.sh
	src/link.sh
	@echo End of compilation
	@echo
	@echo "   and then run 'shw'"
	@echo "------------------------------------------------------------------"

clean: 
	rm -rf bin/*.obj bin/*.o bin/*.mod  
	rm -rf bin/imodel*
	rm -rf *~