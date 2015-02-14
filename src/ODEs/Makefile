
include ../../Makefile.in

CC=gcc
FF=gfortran-4.4
SRC=src/

all:
	$(CC) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O3 -DNDEBUG  $(SRC)"dop853Mex.c"
	$(CC) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O3 -DNDEBUG  $(SRC)"dopri5Mex.c"
	$(CC) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O3 -DNDEBUG  $(SRC)"options.c"
	$(CC) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O3 -DNDEBUG  $(SRC)"tif.c"
	$(FF) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -fexceptions -fbackslash -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 -O3  $(SRC)"dop853.f"
	$(FF) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -fexceptions -fbackslash -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 -O3  $(SRC)"dopri5.f"
	$(FF) -c  -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -fexceptions -fbackslash -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 -O3  $(SRC)"tif_test.f"
	$(CC) -O3 -pthread -shared -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "dop853Mex.mexa64"   "dop853Mex.o"  "options.o"  "tif.o"  "dop853.o"  "tif_test.o"  -lgfortran -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	$(CC) -O3 -pthread -shared -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "dopri5Mex.mexa64"   "dopri5Mex.o"  "options.o"  "tif.o"  "dopri5.o"  "tif_test.o"  -lgfortran -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	rm *.o
clean:
	rm *.mexa64
