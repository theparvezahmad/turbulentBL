# Author		: Parvez Ahmad <pahmed333@gmail.com>
#=============================================================================
 
FC 		= ifort
MPIFC 	= mpiifort
PGRM 	= dnssetbl
RES 	= combineFile
SUF 	= .f90
EXE 	= .x
OBJ 	= .o

FLAGS = -O3 -align all -ansi_alias -fminshared -fno-alias -unroll\
		 -vec_report0 -no-prec-div -warn none -check none
#FLAGS = -traceback -stand f90 -check all -warn all

.PHONY: begin restart cmd1 cmd2

begin:$(PGRM) cmd2

$(PGRM):$(PGRM)$(OBJ)
	$(MPIFC) $(FLAGS) $(PGRM)$(OBJ) -o $(PGRM)$(EXE)

$(PGRM)$(OBJ):$(PGRM)$(SUF)
	$(MPIFC) $(FLAGS) -c $(PGRM)$(SUF)

restart:$(PGRM) $(RES) cmd1 cmd2 # 

$(RES):$(RES)$(OBJ)
	$(FC) $(RES)$(OBJ) -o $(RES)$(EXE)

$(RES)$(OBJ):$(RES)$(SUF)
	$(FC) -c $(RES)$(SUF)

cmd1:
	./$(RES)$(EXE)

cmd2:
	mpdboot
	export I_MPI_PIN_PROCESSOR_LIST=16-23
	nohup mpiexec -np 8 ./$(PGRM)$(EXE) &  

clobber:
	rm -f Output/*
	rm -f Input/3D.dat Input/statData.dat
	rm -f *$(EXE) *$(OBJ) *.out
