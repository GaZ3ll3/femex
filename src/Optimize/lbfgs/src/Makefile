include ../../../../Makefile.in

CXX = g++
FC = gfortran-4.4
OPT = -O3
FFLAGS = -c -I$(MATLAB_ROOT)extern/include\
 -I$(MATLAB_ROOT)simulink/include -fexceptions\
  -fbackslash -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 $(OPT)

F_SRC = $(wildcard *.f)
F_OBJ = $(patsubst %.f, %.o, $(F_SRC))

C_SRC = $(wildcard *.cpp)
C_OBJ = $(patsubst %.cpp, %.o, $(C_SRC))
C_FLAGS =  -c  -I$(MATLAB_ROOT)extern/include \
-I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE\
 -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG

TARGET = lbfgsb.mexa64

LD = $(OPT) -pthread -shared -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o
LINK =  -lgfortran -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm

all: c f 
	$(CXX) $(LD) $(TARGET) $(C_OBJ) $(F_OBJ) $(LINK)
	rm -rf *.o
	mv $(TARGET) ../


f: $(F_SRC)
	$(foreach var, $(F_SRC), $(FC) $(FFLAGS) $(var);)


c: $(C_SRC)
	$(foreach var, $(C_SRC),$(CXX) $(C_FLAGS) $(var);)
	

distclean:
	rm -rf *.o *.mexa64

	