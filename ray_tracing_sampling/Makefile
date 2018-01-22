##
## Receiver function inversion makefile 
## Jan Dettmer jand@uvic.ca 2012
##
## Mac OSX
#FC      = mpif90 -C -L/usr/lib -llapack -lblas
#FC      = mpif90 -fc=pgfortran -L/usr/lib -llapack -lblas -fast -O4 -fastsse
include ./Makefile.compiler

##############################################################################
##
##  Define binary and object paths:
##
BIN = ./bin
OBJ = ./obj
OBJS = $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/alloc_obj.o $(OBJ)/read_input.o $(OBJ)/loglhood.o \
       $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o $(OBJ)/mod_sac_io.o \
       $(OBJ)/four1.o $(OBJ)/fourrow.o $(OBJ)/realft.o $(OBJ)/convlv.o $(OBJ)/quicksort.o \
       $(OBJ)/pythag.o $(OBJ)/prjmh_temper_rf.o
OBJSR= $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/alloc_obj.o $(OBJ)/read_input.o $(OBJ)/loglhood.o \
       $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o $(OBJ)/mod_sac_io.o \
       $(OBJ)/four1.o $(OBJ)/fourrow.o $(OBJ)/realft.o $(OBJ)/convlv.o $(OBJ)/quicksort.o \
       $(OBJ)/pythag.o $(OBJ)/replica.o
FO = ""
LIB =
##############################################################################

all: bin/prjmh_temper_rf

prjmh_temper_rf: bin/prjmh_temper_rf

replica: bin/replica

clean: clean.done
	rm -f *.done
	(cd ./swd/;/bin/rm -f *.o *.mod *.l core *.a;)
clean.done:
	rm -f bin/*
	rm -f *.o *.mod $(OBJ)/*

##############################################################################
#
# Targets for linking of programs:
#
# PRJMH temper Receiver:
#
$(BIN)/prjmh_temper_rf: $(OBJS) $(FO)
	$(FC)  -o $(BIN)/prjmh_temper_rf $(OBJS) $(FO) $(LIB)

$(BIN)/replica: $(OBJSR) $(FO)
	$(FC)  -o $(BIN)/replica $(OBJSR) $(FO) $(LIB)

$(OBJ)/prjmh_temper_rf.o: prjmh_temper_rf.f90 
	$(FC) -c -o $(OBJ)/prjmh_temper_rf.o prjmh_temper_rf.f90 $(LIB)

$(OBJ)/replica.o: replica.f90 
	$(FC) -c -o $(OBJ)/replica.o replica.f90 $(LIB)

$(OBJ)/rjmcmc_com.o: rjmcmc_com.f90 $(OBJ)/data_type.o; $(FC) -c rjmcmc_com.f90 -o $(OBJ)/rjmcmc_com.o $(OBJ)/data_type.o $(LIB)

$(OBJ)/data_type.o: data_type.f90; $(FC) -c data_type.f90 -o $(OBJ)/data_type.o $(LIB)

##############################################################################
#
# Targets for separate subroutines:
##!libswd:
	#(cd ./swd/;make;)
	#/bin/ln -sf ./swd/*.mod .

##############################################################################
#
# Targets to compile NR subroutines that are in separate files:
#
$(OBJ)/read_input.o: read_input.f90
	$(FC) -c -o $(OBJ)/read_input.o read_input.f90

$(OBJ)/alloc_obj.o: alloc_obj.f90
	$(FC) -c -o $(OBJ)/alloc_obj.o alloc_obj.f90

$(OBJ)/loglhood.o: $(OBJ)/quicksort.o loglhood.f90
	$(FC) -c -o $(OBJ)/loglhood.o loglhood.f90

$(OBJ)/ludcmp.o: ludcmp.f90
	$(FC) -c -o $(OBJ)/ludcmp.o ludcmp.f90

$(OBJ)/svdcmp.o: svdcmp.f90
	$(FC) -c -o $(OBJ)/svdcmp.o svdcmp.f90

$(OBJ)/convlv.o: convlv.f90
	$(FC) -c -o $(OBJ)/convlv.o convlv.f90

$(OBJ)/realft.o: realft.f90
	$(FC) -c -o $(OBJ)/realft.o realft.f90

$(OBJ)/four1.o: four1.f90
	$(FC) -c -o $(OBJ)/four1.o four1.f90

$(OBJ)/fourrow.o: fourrow.f90
	$(FC) -c -o $(OBJ)/fourrow.o fourrow.f90

$(OBJ)/pythag.o: pythag.f90
	$(FC) -c -o $(OBJ)/pythag.o pythag.f90

##############################################################################
#
#  Targets for numerical recipes modules:
#
$(OBJ)/nrtype.o: nrtype.f90 
	$(FC) -c -o $(OBJ)/nrtype.o nrtype.f90

$(OBJ)/nr.o: nr.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nr.o nr.f90

$(OBJ)/nrutil.o: nrutil.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nrutil.o nrutil.f90

$(OBJ)/quicksort.o: quicksort.f90
	$(FC) -c -o $(OBJ)/quicksort.o quicksort.f90

$(OBJ)/mod_sac_io.o : mod_sac_io.f90
	$(FC) -c -o $(OBJ)/mod_sac_io.o mod_sac_io.f90

$(OBJ)/subroutineR-quiet.o: $(OBJ)/raytracing/subroutineR-quiet.f90
	$(FC) -c -o $(OBJ)/raytracing/subroutineR-quiet.o ./raytracing/subroutineR-quiet.f90

##############################################################################
# ...this is the end my fiend
#
# EOF
