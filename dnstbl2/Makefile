############################ -*- Mode: Makefile -*- ###########################
## Usage make; make ;make clean
## Author          : <Fuaad P A  fukipa@gmail.com>
###############################################################################
# No modifications after this line.

FFLAGS_OPT	=  -O3 -align all -ansi_alias -fminshared -warn usage -check none\
              -fno-alias -unroll -vec_report0 -no-prec-div 
#FFLAGS_OPT	= -O0 -g -traceback -fpe:0 -check all 

MPIFC      = mpiifort
FC      = ifort
FFLAGS  = -O3
EXES    = *.x  



SRCMOD	=dnstbl.f90
AVGMOD	=averages.f90
OBJMOD	=dnstbl.o

.SUFFIXES: .f90 .o

all: 
	$(MPIFC) $(FFLAGS_OPT) $(SRCMOD) -o dnstbl.x 

averages: 
	$(FC) $(FFLAGS) $(AVGMOD) -o averages.x  
clean:
	rm -rf *.o *.mod
clobber: clean
	rm -f *.out *.tp means medians taustress *.xy
