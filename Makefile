# -----------------------------------------------------------------
# Version: 1.0
# Date: Nov 2019
# Makefile for SHUD v 1.0
# -----------------------------------------------------------------
# Programmer: Lele Shu (lele.shu@gmail.com)
# SHUD model is a heritage of Penn State Integrated Hydrologic Model (PIHM).
# -----------------------------------------------------------------
#  Prerequisite:
#  1 install sundials 5.0+ via https://computation.llnl.gov/projects/sundials/sundials-software.
#  2 If parallel-computing is prefered, please install OpenMP.
#	 For mac: 
#	  		brew install llvm clang
#			brew install libomp
#			compile flags for OpenMP: 
#				-Xpreprocessor -fopenmp -lomp
#			Library/Include paths:
#				-L/usr/local/opt/libomp/lib 
#				-I/usr/local/opt/libomp/include
#			
# -----------------------------------------------------------------
# Configure this File:
# 1 Path of SUNDIALS_DIR. [CRITICAL]
# 2 Path of OpenMP if parallel is preffered.
# 3 Path of SRC_DIR, default is "SRC_DIR = ."
# 4 Path of BUILT_DIR, default is "BUILT_DIR = ."
# -----------------------------------------------------------------
SUNDIALS_DIR = ~/sundials
# SUNDIALS_DIR = /usr/local/sundials


SHELL = /bin/sh
BUILDDIR = .
SRC_DIR = src

LIB_USR = ${SUNDIALS_DIR}/lib
LIB_SYS = /usr/local/lib/
INC_OMP = /usr/local/opt/libomp/include
LIB_OMP = /usr/local/opt/libomp/lib

INC_MPI = /usr/local/opt/open-mpi

TARGET_EXEC     = ${BUILDDIR}/shud
TARGET_OMP      = ${BUILDDIR}/shud_omp
TARGET_DEBUG    = ${BUILDDIR}/shud_debug

MAIN_shud 		= ${SRC_DIR}/SHUD.cpp
MAIN_OMP 		= ${SRC_DIR}/SHUD.cpp
MAIN_DEBUG 		= ${SRC_DIR}/SHUD.cpp

# If compile on Cluster
# CC       = g++
# MPICC    = mpic++
# LK_OMP   = -fopenmp -lsundials_nvecopenmp
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SUNDIALS_DIR}/lib

CC       = /usr/bin/g++
MPICC    = /usr/local/bin/mpic++
CFLAGS   = -O3 -g  -std=c++11
#STCFLAG     = -static
LDFLAGS  = -Wl, -rpath ${SUNDIALS_DIR}/lib
LIBS     = -lm
SRC    	= ${SRC_DIR}/classes/*.cpp \
		  ${SRC_DIR}/ModelData/*.cpp \
		  ${SRC_DIR}/Model/*.cpp \
		  ${SRC_DIR}/Equations/*.cpp

SRC_H	= ${SRC_DIR}/classes/*.hpp \
		  ${SRC_DIR}/ModelData/*.hpp \
		  ${SRC_DIR}/Model/*.hpp \
		  ${SRC_DIR}/Equations/*.hpp


INCLUDES = -I ${SUNDIALS_DIR}/include \
		   -I ${INC_OMP} \
		   -I ${SRC_DIR}/Model \
		   -I ${SRC_DIR}/ModelData \
		   -I ${SRC_DIR}/classes \
		   -I ${SRC_DIR}/Equations 

		  
LIBRARIES = -L ${LIB_OMP} \
			-L ${SUNDIALS_DIR}/lib \
			-L ${LIB_SYS}

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial
LK_OMP	= -Xpreprocessor -fopenmp -lomp -lsundials_nvecopenmp



all:
	make clean
	make shud
	@echo
check:
	./shud
	@echo
help:
	@(echo)
	@echo "Usage:"
	@(echo '       make all	    	- make both shud and shud_omp')
	@(echo '       make cvode	    - install SUNDIALS/CVODE to ~/sundials')
	@(echo '       make shud     	- make shud executable')
	@(echo '       make shud_omp    - make shud_omp with OpenMP support')
	@(echo)
	@(echo '       make clean    	- remove all executable files')
	@(echo)
cvode CVODE:
	@echo '...Install SUNDIALS/CVODE for your ...'
	chmod +x ./script/installSundials.sh
	./script/installSundials.sh
	
	@echo 

shud SHUD: ${MAIN_shud} $(SRC) $(SRC_H)
	@echo '...Compiling shud ...'
	@echo $(CC) $(CFLAGS) ${STCFLAG} ${INCLUDES} ${LIBRARIES} ${LDFLAGS} -o ${TARGET_EXEC} ${MAIN_shud} $(SRC)  $(LK_FLAGS)
	@echo
	@echo
	$(CC) $(CFLAGS) ${INCLUDES} ${STCFLAG} ${LIBRARIES} ${LDFLAGS} -o ${TARGET_EXEC} ${MAIN_shud} $(SRC)  $(LK_FLAGS)
	@echo
	@echo
	@echo " ${TARGET_EXEC} is compiled successfully!"
	@echo

shud_omp: ${MAIN_OMP}  $(SRC) $(SRC_H)
	@echo '...Compiling shud_OpenMP ...'
	@echo $(CC) $(CFLAGS) ${STCFLAG} ${LDFLAGS} -D_OPENMP_ON ${INCLUDES} ${LIBRARIES} -o ${TARGET_OMP}   ${MAIN_OMP} $(SRC)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo
	$(CC) $(CFLAGS)  ${STCFLAG} ${LDFLAGS} -D_OPENMP_ON ${INCLUDES} ${LIBRARIES} -o ${TARGET_OMP}   ${MAIN_OMP} $(SRC)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo " ${TARGET_OMP} is compiled successfully!"
	@echo
	@echo

clean:
	@echo "Cleaning ... "
	@echo
	@echo "  rm -f *.o"
	@rm -f *.o
	
	@echo "  rm -f ${TARGET_EXEC}"
	@rm -f ${TARGET_EXEC}
	
	@echo "  rm -f ${TARGET_OMP}"
	@rm -f ${TARGET_OMP}
	
	@echo
	@echo "Done."
	@echo





