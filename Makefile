CXX = g++-4.8
CC  = gcc-4.8
FF  = gfortran-4.8
Opt = -Ofast

include Makefile.in

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			-DMX_COMPAT_32 $(Opt) -DNDEBUG -fopenmp -ffast-math

CXX_INCLUDE = -I./include/mexplus \
			  -I./include/pprint \
			  -I./include/Utility \
			  -I./include/exprtk \
			  -I$(MATLAB_ROOT)extern/include \
			  -I$(MATLAB_ROOT)simulink/include


MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map \
			   -Wl,--no-undefined 
			   
CXX_LIBS = -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 \
		   -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp
		   
	
	
# Common use
########################################################

	
########################################################		   
		   
SRC = src/

# mesh  build
MESH = $(SRC)Mesh/private/
TRIANGLELIB = $(MESH)triangle/
MESH_SRCS = $(wildcard $(MESH)*.cc)
MESH_OBJS = $(patsubst $(MESH)%.cc, %.mesh.o, $(MESH_SRCS))
MESH_BINS = $(patsubst $(MESH)%.cc, $(MESH)%_.mexa64, $(MESH_SRCS))

%.mesh.o: $(MESH)%.cc 
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(MESH)%_.mexa64: %.mesh.o $(TRIANGLELIB)triangle.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(TRIANGLELIB)triangle.o $(CXX_LIBS)
	
	
# Triangle build

CSWITCHES = $(Opt) -DLINUX -fPIC -shared -I/usr/X11R6/include -L/usr/X11R6/lib
TRILIBDEFS = -DTRILIBRARY

$(TRIANGLELIB)triangle.o:$(TRIANGLELIB)triangle.c $(TRIANGLELIB)triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c -o $(TRIANGLELIB)triangle.o \
	$(TRIANGLELIB)triangle.c 

# Integrator build

INT = $(SRC)Integrator/private/
INT_SRCS = $(wildcard $(INT)*.cc)
INT_OBJS = $(patsubst $(INT)%.cc, %.int.o, $(INT_SRCS))
INT_BINS = $(patsubst $(INT)%.cc, $(INT)%_.mexa64, $(INT_SRCS))

%.int.o: $(INT)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(INT)%_.mexa64: %.int.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)


# Assembler build

ASR = $(SRC)Assembler/private/
ASR_SRCS = $(wildcard $(ASR)*.cc)
ASR_OBJS = $(patsubst $(ASR)%.cc, %.asr.o, $(ASR_SRCS))
ASR_BINS = $(patsubst $(ASR)%.cc, $(ASR)%_.mexa64, $(ASR_SRCS))

%.asr.o: $(ASR)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
$(ASR)%_.mexa64: %.asr.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)



BOD = $(SRC)Boundary/private/
BOD_SRCS = $(wildcard $(BOD)*.cc)
BOD_OBJS = $(patsubst $(BOD)%.cc, %.bod.o, $(BOD_SRCS))
BOD_BINS = $(patsubst $(BOD)%.cc, $(BOD)%_.mexa64, $(BOD_SRCS))

%.bod.o: $(BOD)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(BOD)%_.mexa64: %.bod.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)
	
	
SLR = $(SRC)Solver/private/
SLR_SRCS = $(wildcard $(SLR)*.cc)
SLR_OBJS = $(patsubst $(SLR)%.cc, %.slr.o, $(SLR_SRCS))
SLR_BINS = $(patsubst $(SLR)%.cc, $(SLR)%_.mexa64, $(SLR_SRCS))

%.slr.o: $(SLR)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
$(SLR)%_.mexa64: %.slr.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)
	
	
	
	
SVR = $(SRC)PDEcs/private/
SVR_SRCS = $(wildcard $(SVR)*.cc)
SVR_OBJS = $(patsubst $(SVR)%.cc, %.svr.o, $(SVR_SRCS))
SVR_BINS = $(patsubst $(SVR)%.cc, $(SVR)%_.mexa64, $(SVR_SRCS))

%.svr.o: $(SVR)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
$(SVR)%_.mexa64: %.svr.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)
##############################################################
# modules
##########

MOD = modules/

ASE= $(MOD)AssembleExtension/private/
ASE_SRCS = $(wildcard $(ASE)*.cc)
ASE_OBJS = $(patsubst $(ASE)%.cc, %.ase.o, $(ASE_SRCS))
ASE_BINS = $(patsubst $(ASE)%.cc, $(ASE)%_.mexa64, $(ASE_SRCS))

%.ase.o: $(ASE)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(ASE)%_.mexa64: %.ase.o 
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)

###########
MSE= $(MOD)MeshExtension/private/
MSE_SRCS = $(wildcard $(MSE)*.cc)
MSE_OBJS = $(patsubst $(MSE)%.cc, %.mse.o, $(MSE_SRCS))
MSE_BINS = $(patsubst $(MSE)%.cc, $(MSE)%_.mexa64, $(MSE_SRCS))

%.mse.o: $(MSE)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(MSE)%_.mexa64: %.mse.o 
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)

##############################################################
# ILUPACK make
ILUPACK_ROOT = ./$(SRC)Solver/ilupack/
ILUPACK = ./$(SRC)Solver/ilupack/matlabsrc/
ILUPACK_PATH = ./$(SRC)Solver/ilupack/mex/
ILUPACK_FLAGS = -c  -I$(ILUPACK_ROOT)include -I$(MATLAB_ROOT)extern/include -I$(MATLAB_ROOT)simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -D_LONG_INTEGER_ -D_MUMPS_MATCHING_ -D__UNDERSCORE__ $(Opt) -DNDEBUG  

ILUPACK_LINKS =  $(Opt) -pthread -shared -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o 
ILUPACK_LIBS  =  -L$(ILUPACK_ROOT)lib/GNU64_long -lilupack -lmumps -lamd -lsparspak -lblaslike -lmwlapack -lmwblas -lmetis -lm -lc -lgfortran -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++


ILUPACK_SRCS = $(wildcard $(ILUPACK)*.c)
ILUPACK_OBJS = $(patsubst $(ILUPACK)%.c, $(ILUPACK_PATH)%.o, $(ILUPACK_SRCS))
ILUPACK_BINS = $(patsubst $(ILUPACK)%.c, $(ILUPACK_PATH)%.mexa64, $(ILUPACK_SRCS))

$(ILUPACK_PATH)%.o: $(ILUPACK)%.c
	$(CC) $(ILUPACK_FLAGS) $< -o $@
$(ILUPACK_PATH)%.mexa64: $(ILUPACK_PATH)%.o
	$(FF) $(ILUPACK_LINKS) $@ $< $(ILUPACK_LIBS)

##############################################################	
# The action starts here.
all: $(MESH_BINS) $(ASR_BINS) $(INT_BINS) $(BOD_BINS) $(SLR_BINS) $(SVR_BINS) \
     $(ASE_BINS) $(MSE_BINS) $(ILUPACK_BINS)

module: $(ASE_BINS) $(MSE_BINS)

distclean:
	rm -rf $(MESH)*_.mexa64 $(INT)*_.mexa64 $(ASR)*_.mexa64 $(BOD)*_.mexa64 $(SLR)*_.mexa64 $(ASE)*_.mexa64 $(MSE)*_.mexa64 $(TRIANGLELIB)triangle.o $(ILUPACK_PATH)*.mexa64 

clean:
	rm -rf $(MESH)*_.mexa64 $(INT)*_.mexa64 $(ASR)*_.mexa64 $(BOD)*_.mexa64 $(SLR)*_.mexa64 $(ASE)*_.mexa64 $(MSE)*_.mexa64 $(TRIANGLELIB)triangle.o
