CXX = g++
CC  = gcc
FF  = gfortran
Opt = -Ofast

include Makefile.in

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			-DMX_COMPAT_32 $(Opt) -DNDEBUG -fopenmp -ffast-math 

CXX_INCLUDE = -I./include/mexplus \
			  -I./include/pprint \
			  -I./include/Utility \
			  -I./include/exprtk \
                          -I/usr/include/eigen3\
			  -I$(MATLAB_ROOT)extern/include \
			  -I$(MATLAB_ROOT)simulink/include


MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map \
			   -Wl,--no-undefined -lblas -llapack
			   
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

# Mesh3
MESH3 = $(SRC)Mesh3/private/
TETGENLIB = $(MESH3)tetgen/
MESH3_SRCS = $(wildcard $(MESH3)*.cc)
MESH3_OBJS = $(patsubst $(MESH3)%.cc, %.mesh3.o, $(MESH3_SRCS))
MESH3_BINS = $(patsubst $(MESH3)%.cc, $(MESH3)%_.mexa64, $(MESH3_SRCS))

%.mesh3.o: $(MESH3)%.cc 
	$(CXX) -c $(TETLIBDEFS) $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(MESH3)%_.mexa64: %.mesh3.o $(TETGENLIB)tetgen.o $(TETGENLIB)predicates.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(TETGENLIB)tetgen.o $(TETGENLIB)predicates.o $(CXX_LIBS)

# Tetgen build

CSWITHCES3 = $(Opt) -DLINUX -fPIC -shared 
TETLIBDEFS = -DTETLIBRARY

$(TETGENLIB)tetgen.o:$(TETGENLIB)tetgen.cxx
	$(CXX) $(CSWITHCES3) $(TETLIBDEFS) -c -o $(TETGENLIB)tetgen.o \
	$(TETGENLIB)tetgen.cxx 
	
$(TETGENLIB)predicates.o:$(TETGENLIB)predicates.cxx
	$(CXX) $(CSWITHCES3) $(TETLIBDEFS) -c -o $(TETGENLIB)predicates.o \
	$(TETGENLIB)predicates.cxx 

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
	
	
DOM = $(SRC)DOM/private/
DOM_SRCS = $(wildcard $(DOM)*.cc)
DOM_OBJS = $(patsubst $(DOM)%.cc, %.dom.o, $(DOM_SRCS))
DOM_BINS = $(patsubst $(DOM)%.cc, $(DOM)%_.mexa64, $(DOM_SRCS))

%.dom.o: $(DOM)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
$(DOM)%_.mexa64: %.dom.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)

CEL = $(SRC)FMM/private/
CEL_SRCS = $(wildcard $(CEL)*.cc)
CEL_OBJS = $(patsubst $(CEL)%.cc, %.cel.o, $(CEL_SRCS))
CEL_BINS = $(patsubst $(CEL)%.cc, $(CEL)%_.mexa64, $(CEL_SRCS))

%.cel.o: $(CEL)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
$(CEL)%_.mexa64: %.cel.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)
	
	
QUA = $(SRC)/QuadTree/private/
QUA_BINS = $(QUA)QuadTree_.mexa64

$(QUA)QuadTree_.mexa64: $(QUA)Quadtree_.cc $(QUA)QuadTree.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(QUA)Quadtree_.cc -o Quadtree_.qua.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(QUA)QuadTree.cc -o QuadTree.qua.o
	$(CXX) $(MATLAB_LINKS) -o $@ Quadtree_.qua.o QuadTree.qua.o $(CXX_LIBS)
	
TRE = $(SRC)/Treecode/private/
TRE_BINS = $(TRE)Treecode_.mexa64
	
$(TRE)Treecode_.mexa64: $(TRE)quadtree.cpp $(TRE)treecode.cpp $(TRE)Treecode_.cpp $(TRE)treecode.h
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(TRE)quadtree.cpp -o quadtree.tre.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(TRE)treecode.cpp -o treecode.tre.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(TRE)Treecode_.cpp -o Treecode_.tre.o
	$(CXX) $(MATLAB_LINKS) -o $@ Treecode_.tre.o treecode.tre.o quadtree.tre.o $(CXX_LIBS)
	
	
	

RAD = $(SRC)Radfmm/private/
RAD_SRCS = $(wildcard $(RAD)*.cpp)
RAD_OBJS = $(patsubst $(RAD)%.cpp, %.rad.o, $(RAD_SRCS))
RADK_BINS = $(RAD)Radfmmk_.mexa64


$(RAD)Radfmmk_.mexa64: $(RAD)H2_2D_Node.cpp $(RAD)H2_2D_Tree.cpp $(RAD)kernel_Base.cpp $(RAD)kernel_Types.cpp $(RAD)Radfmmk.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(RAD)H2_2D_Node.cpp -o H2_2D_Node.rad.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(RAD)H2_2D_Tree.cpp -o H2_2D_Tree.rad.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(RAD)kernel_Base.cpp -o kernel_Base.rad.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(RAD)kernel_Types.cpp -o kernel_Types.rad.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(RAD)Radfmmk.cpp -o Radfmmk.rad.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.rad.o H2_2D_Tree.rad.o kernel_Base.rad.o kernel_Types.rad.o  Radfmmk.rad.o $(CXX_LIBS)
	
	

KIF = $(SRC)KIFMM2D/private/
KIF_SRCS = $(wildcard $(KIF)*.cpp)
KIF_OBJS = $(patsubst $(KIF)%.cpp, %.kif.o, $(KIF_SRCS))
KIF_BINS = $(KIF)kifmm_.mexa64

$(KIF)kifmm_.mexa64: $(KIF)LET.cpp $(KIF)LET_Node.cpp $(KIF)kifmm.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(KIF)LET_Node.cpp -o LET_Node.kif.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(KIF)LET.cpp -o LET.kif.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(KIF)kifmm.cpp -o kifmm.kif.o
	$(CXX) $(MATLAB_LINKS) -o $@ LET_Node.kif.o LET.kif.o kifmm.kif.o $(CXX_LIBS)

EDN = $(SRC)Eddington/private/
EDN_SRCS=$(wildcard $(EDN)*.cpp)
EDN_OBJS=$(patsubst $(EDN)%.cpp, %.edn.o, $(EDN_SRCS))
EDN_UU_BINS=$(EDN)RadfmmUU_.mexa64
EDN_UC_BINS=$(EDN)RadfmmUC_.mexa64
EDN_US_BINS=$(EDN)RadfmmUS_.mexa64
EDN_CC_BINS=$(EDN)RadfmmCC_.mexa64
EDN_CS_BINS=$(EDN)RadfmmCS_.mexa64
EDN_SS_BINS=$(EDN)RadfmmSS_.mexa64	


$(EDN)RadfmmUU_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmUU.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmUU.cpp -o RadfmmUU.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmUU.edn.o $(CXX_LIBS)


$(EDN)RadfmmUC_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmUC.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmUC.cpp -o RadfmmUC.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmUC.edn.o $(CXX_LIBS)
	
	
$(EDN)RadfmmUS_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmUS.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmUS.cpp -o RadfmmUS.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmUS.edn.o $(CXX_LIBS)
	

$(EDN)RadfmmCC_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmCC.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmCC.cpp -o RadfmmCC.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmCC.edn.o $(CXX_LIBS)


$(EDN)RadfmmCS_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmCS.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmCS.cpp -o RadfmmCS.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmCS.edn.o $(CXX_LIBS)


$(EDN)RadfmmSS_.mexa64: $(EDN)H2_2D_Node.cpp $(EDN)H2_2D_Tree.cpp $(EDN)kernel_Base.cpp $(EDN)kernel_Types.cpp $(EDN)RadfmmSS.cpp
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Node.cpp -o H2_2D_Node.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)H2_2D_Tree.cpp -o H2_2D_Tree.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Base.cpp -o kernel_Base.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)kernel_Types.cpp -o kernel_Types.edn.o
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(EDN)RadfmmSS.cpp -o RadfmmSS.edn.o
	$(CXX) $(MATLAB_LINKS) -o $@ H2_2D_Node.edn.o H2_2D_Tree.edn.o kernel_Base.edn.o kernel_Types.edn.o  RadfmmSS.edn.o $(CXX_LIBS)

MCR = $(SRC)/MCRT/private/
MCR_BINS = $(MCR)MCRT_.mexa64
	
$(MCR)MCRT_.mexa64: $(MCR)MCRT.cc 
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $(MCR)MCRT.cc -o MCRT.mcr.o
	$(CXX) $(MATLAB_LINKS) -o $@ MCRT.mcr.o  $(CXX_LIBS)
	
#ADJ = $(SRC)Adjoint/private/
#ADJ_SRCS = $(wildcard $(ADJ)*.cc)
#ADJ_OBJS = $(patsubst $(ADJ)%.cc, %.adj.o, $(ADJ_SRCS))
#ADJ_BINS = $(patsubst $(ADJ)%.cc, $(ADJ)%_.mexa64, $(ADJ_SRCS))

#%.adj.o: $(ADJ)%.cc
#	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
#$(ADJ)%_.mexa64: %.adj.o
#	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS)
##############################################################

Optimize=$(SRC)Optimize/

Ipopt: 
	cd $(Optimize) && make 

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
all: $(MESH_BINS) $(MESH3_BINS) $(ASR_BINS) $(INT_BINS) $(BOD_BINS) $(SLR_BINS) $(ILUPACK_BINS) \
$(DOM_BINS) $(CEL_BINS) $(QUA_BINS) $(TRE_BINS) $(RAD_BINS) $(RADK_BINS) $(MCR_BINS) $(EDN_UU_BINS) \
$(EDN_UC_BINS) $(EDN_US_BINS) $(EDN_CC_BINS) $(EDN_CS_BINS) $(EDN_SS_BINS) $(KIF_BINS)
	rm -rf $(TRIANGLELIB)triangle.o\
	rm -rf *.o

distclean:
	rm -rf $(MESH)*_.mexa64 $(INT)*_.mexa64 $(ASR)*_.mexa64 $(BOD)*_.mexa64 $(SLR)*_.mexa64 $(ASE)*_.mexa64 \
	 $(MSE)*_.mexa64 $(TRIANGLELIB)triangle.o $(ILUPACK_PATH)*.mexa64 $(Optimize)ipopt/ipopt.mexa64 \
	 $(DOM)*_.mexa64 $(CEL)*_.mexa64 $(QUA)*_.mexa64 $(TRE)*_.mexa64 $(RAD)*_.mexa64 $(MCR)*_.mexa64

clean:
	rm -rf $(MESH)*_.mexa64 $(INT)*_.mexa64 $(ASR)*_.mexa64 $(BOD)*_.mexa64 $(SLR)*_.mexa64 \
	$(ASE)*_.mexa64 $(MSE)*_.mexa64  $(DOM)*_.mexa64 $(CEL)*_.mexa64 $(TRIANGLELIB)triangle.o \
	$(QUA)*_.mexa64 $(TRE)*_.mexa64 $(RAD)*_.mexa64 $(MCR)*_.mexa64
