# Only works on 64-bit.
# Edit these variables.
MEXSUFFIX   = mexa64
MEX         = /usr/local/MATLAB/R2014b/bin/mex
CXX         = g++
PARDISOHOME = ../

MATLAB_ROOT = /usr/local/MATLAB/R2014b/

Opt = -Ofast

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
                        -D_GNU_SOURCE -fno-omit-frame-pointer -Wno-write-strings -pthread\
                        -DMX_COMPAT_32 $(Opt) -DNDEBUG -fopenmp -ffast-math

CXX_INCLUDE = -I$(MATLAB_ROOT)extern/include \
                          -I$(MATLAB_ROOT)simulink/include


MATLAB_LINKS = $(Opt) -pthread \
                           -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map \
                           -Wl,--no-undefined -lblas -llapack

CXX_LIBS = -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 \
                   -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp


# Do not edit anything past here (unless you know what you're doing).
#MEXFLAGS   = -cxx CXX=$(CXX) CC=$(CXX) LD=$(CXX) -largeArrayDims \
#             -lmwlapack -lmwblas \
#             -L../libpardiso500-GNU472-X86-64.so -L/usr/lib/x86_64-linux-gnu/libgomp.so.1 -lgfortran -lpthread -lm -L../libpardiso500-GNU472-X86-64.so 

TARGETS_WO_MEXSUFFIX = pardisoinit pardisoreorder \
          pardisofactor pardisosolve \
          pardisofree

TARGETS = pardisoinit.$(MEXSUFFIX) pardisoreorder.$(MEXSUFFIX) \
          pardisofactor.$(MEXSUFFIX) pardisosolve.$(MEXSUFFIX) \
          pardisofree.$(MEXSUFFIX)

all: $(TARGETS_WO_MEXSUFFIX)

%: %.cpp common.cpp matlabmatrix.cpp sparsematrix.cpp \
                pardisoinfo.cpp
	$(CXX) $(CXX_INCLUDE) -o $@ $< $(MATLAB_LINKS) $(CXX_LIBS)

clean:
	rm -f $(TARGETS) 
