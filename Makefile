CXX = g++

Opt = -Ofast

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			-DMX_COMPAT_32 $(Opt) -DNDEBUG -fopenmp

CXX_INCLUDE = -I./include/mexplus \
			  -I./include/pprint \
			  -I./include/Utility \
			  -I/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/include \
			  -I/usr/local/MATLAB/MATLAB_Production_Server/R2013a/simulink/include


MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/lib/glnxa64/mexFunction.map \
			   -Wl,--no-undefined 
			   
CXX_LIBS = -Wl,-rpath-link,/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64 \
		   -L/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp
		   
	
	
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
	
# The action starts here.

# all: $(BIN)MeshGen.mexa64 $(SRC)triangle.o $(BIN)MatrixAssem.mexa64 $(BIN)BCAssem.mexa64 clean

all: $(MESH_BINS) $(ASR_BINS) $(INT_BINS) $(BOD_BINS)

clean:
	rm -rf $(MESH)*_.mexa64 $(INT)*_.mexa64 $(ASR)*_.mexa64 $(BOD)*_.mexa64 $(TRIANGLELIB)triangle.o