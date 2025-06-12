# makefile for run_sigma_ed4 (replace run_sigma_ed4.mak)


OBJS    = help.o dwavefunctionc_pp_pn.o bystricky.o fpn_chex.o f_NN_sd.o fpn.o b_function.o j_e.o j_N_matrix_s_onof.o sigma_ed.o run_sigma_ed4.f90


F90     = gfortran



FFLAGS  = -fno-automatic -O2




#----------------------------------------------------------
all: $(OBJS) run_sigma_ed4

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

%.o: %.f
	$(F90) $(FFLAGS) -c $<


run_sigma_ed4: 
	$(F90)  $(OBJS) -o run_sigma_ed4



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o run_sigma_ed4




