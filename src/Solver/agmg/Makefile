include ../../../Makefile.in

CC = gcc
# thanks to ABI
F4 = gfortran-4.4
F8 = gfortran-4.8

ALL_C_FLAGS = -c -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE\
  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -Ofast -DNDEBUG
# suppress all warnings
ALL_F_FLAGS = -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include  -fexceptions -fbackslash\
 -fPIC -fno-omit-frame-pointer -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -w -O3 \
 -fforce-addr -fstrength-reduce -fcaller-saves
 
PRE_LD = -O3 -pthread -shared -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  
POST_LD = -lmwlapack -lmwblas -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 \
-lmx -lmex -lmat -lm -lstdc++ -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 \
-lmx -lmex -lmat -lm

all: d z w clean
d: dmtlagmg.o   dagmg_mex.o dagmg_mumps_mex.o
	$(F8) $(PRE_LD) "dmtlagmg.mexa64"  dmtlagmg.o dagmg_mex.o dagmg_mumps_mex.o $(POST_LD)
dmtlagmg.o:
	$(CC) $(ALL_C_FLAGS) "dmtlagmg.c"
dagmg_mex.o:
	$(F4) $(ALL_F_FLAGS) "dagmg_mex.f90"
dagmg_mumps_mex.o:
	$(F4) $(ALL_F_FLAGS) "dagmg_mumps_mex.f90"


	 
z: zmtlagmg.o zagmg_mex.o zagmg_mumps_mex.o
	$(F8) $(PRE_LD) "zmtlagmg.mexa64"  zmtlagmg.o zagmg_mex.o zagmg_mumps_mex.o $(POST_LD)
zmtlagmg.o:
	$(CC) $(ALL_C_FLAGS) "zmtlagmg.c"
zagmg_mex.o:
	$(F4) $(ALL_F_FLAGS) "zagmg_mex.f90"
zagmg_mumps_mex.o:
	$(F4) $(ALL_F_FLAGS) "zagmg_mumps_mex.f90"


w:  dmtlagtwolev.o dagtwolev_mex.o
	$(F8) $(PRE_LD) "dmtlagtwolev.mexa64"  dmtlagtwolev.o dagtwolev_mex.o $(POST_LD)
dmtlagtwolev.o:
	$(CC) $(ALL_C_FLAGS) "dmtlagtwolev.c"
dagtwolev_mex.o:
	$(F4) $(ALL_F_FLAGS) "dagtwolev_mex.f90"

	
	
distclean:
	rm -rf *.o *.mod *.mexa64
	
clean:
	rm -rf *.o *.mod
