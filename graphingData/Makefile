#####################################
# Makefile for GraphPACK
#
# Ryan A. Rossi, Nesreen Ahmed
# Copyright 2012-2017
#####################################

.KEEP_STATE:

ROOTDIR  = $(shell pwd)
OPSYS    = $(shell uname -s )

OPENMP_LIB ?= 1

all: graphpack
	
#CFLAGS = -Wformat -O3 
CFLAGS = -O3 
ifeq ($(OPENMP_LIB),1)
        CFLAGS += -fopenmp
else
        CFLAGS += -Wno-unknown-pragmas
endif
#CFLAGS += -D_GLIBCXX_PARALLEL

# UPDATED 16/11/15
#CFLAGS += -std=c++11
#CFLAGS += -m64 
CFLAGS += -fpeel-loops -ftracer -funit-at-a-time 
#CFLAGS += -floop-parallelize-all 
#-ftree-loop-distribution
CFLAGS += -ffast-math -funroll-loops -fno-strict-aliasing \
	 -fomit-frame-pointer -fexpensive-optimizations -funroll-loops \
	 -fmove-loop-invariants -fprefetch-loop-arrays -ftree-loop-optimize \
	 -ftree-vect-loop-version -ftree-vectorize
	 
CXX          = g++
H_FILES     = graphpack.h

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

IO_SRC 				   = graphpack_utils.cpp graphpack_graph.cpp graphpack_clique_utils.cpp 
MCPACK_SRC 			   = graphpack_heu.cpp graphpack_mc.cpp 
BOUND_LIB_SRC 		   = bound_cores.cpp bound_triangles.cpp bound_triangle_cores.cpp 

MCPACK_MAIN			   = graphpack_driver.cpp


OBJ_MCPACK	= $(MCPACK_MAIN:%.cpp=%.o) $(IO_SRC) $(MCPACK_SRC) $(BOUND_LIB_SRC)
$(OBJ_MCPACK): $(H_FILES) Makefile
graphpack: $(OBJ_MCPACK) $(H_FILES)
	$(CXX) $(CFLAGS) -o graphpack $(OBJ_MCPACK)
		
help:  
	@echo " "
	@echo "Operating System: $(OPSYS) "
	@echo " "
	@echo "USAGE: "
	@echo "make help       To get this listing"
	@echo "make            To compile mcpack"
	@echo "make clean      Remove *.o and executable files"
	@echo "make list       List the compilers in current environment"
	@echo " "

list:
	@echo
	@echo "OPSYS:       $(OPSYS)"
	@echo "ROOTDIR:     $(ROOTDIR)"
	@echo " "

docs:
	@echo "creating documentation in ./html directory"
	doxygen Doxyfile
	
tar: 
	make clean; cd ../;  tar cvzf graphpack.tgz graphpack/
	
clean:
	@echo "Removing *.o and executable files"
	rm -rf *.o graphpack