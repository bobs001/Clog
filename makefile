.SUFFIXES: .c .f .f90 .o .s .mod

NAME = dedx

F90 = gfortran
F77 = gfortran
CC = gcc

#F90FLAGS = -O3  -g -fdefault-real-8
#F90FLAGS = -g -fdefault-real-8
#F90FLAGS = -O3 -g3 -fdefault-real-8 -Wall -ffast-math
F90FLAGS = -g3 -fdefault-real-8 -Wall -ffast-math 

LINK =
#F90FLAGS = -O3  -g -r8 -pg
#LINK = -pg

F77FLAGS = -g -fdefault-real-8
CFLAGS = 

#OBJ_FILES =  ferf.o vars.o mathvars.o physvars.o common.o rate.o acoeff.o dedx.o dedx.models.o main.range.o
OBJ_FILES =  ferf.o vars.o mathvars.o physvars.o common.o rate.o acoeff.o dedx.o dedx.models.o main.acoeff.dedx.o
#OBJ_FILES =  ferf.o vars.o mathvars.o physvars.o common.o rate.o acoeff.o dedx.o dedx.models.o main.dedx.o
#OBJ_FILES =  ferf.o vars.o mathvars.o physvars.o common.o rate.o acoeff.o dedx.o dedx.models.o main.table.o
#OBJ_FILES =  ferf.o vars.o mathvars.o physvars.o common.o rate.o acoeff.o dedx.o dedx.models.o main.range.o


#   1. main.dedx.f90          <== template for stopping power
#      plot_dedx.py           <== python plotting file
#   2. main_acoeff.f90        <== computes A-coeffs and asymptotic limits
#      main.acoeff.dat        <== data file generated by main_acoeff.f90
#      main.acoeff.highE.dat  <== data file generated by main_acoeff.f90
#      main.acoeff.smallE.dat <== data file generated by main_acoeff.f90
#   3. main_acoeff_dedx.f90   <== computes dedx directly and from the A-coeffs
#   4. main_CeI.f90           <== computes C_{eI} [general and Born]
#
#   5. main_symmetry.f90      <== still working on this one: checks symmetry
#                                 of C_{ab} and C_{eI}




INC = mathvars.f90 physvars.f90 vars.f90 

.f90.o: ; $(F90) $(F90FLAGS) -c $*.f90
.f.o:   ; $(F77) $(F77FLAGS) -c $*.f
.c.o:   ; $(CC)  $(CFLAGS)   -c $*.c

$(NAME): $(OBJ_FILES)
	$(F90) $(OBJ_FILES) $(FFLAGS) $(LINK)   -o $(NAME)

$(OBJ_FILES): $(INC) makefile



clean:
	-rm *.o *.mod

clobber:
	-rm *.o  *.mod $(NAME)
