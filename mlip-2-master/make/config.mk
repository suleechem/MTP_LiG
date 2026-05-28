# Directories
PREFIX = .
BIN_DIR = $(PREFIX)/bin
LIB_DIR = $(PREFIX)/lib
FORTRAN_DIR = /usr/local/gfortran/lib

# Compilers for the executable
CC_EXE  = gcc
CXX_EXE = g++
FC_EXE  = gfortran

# Compilers for the library
CC_LIB  = gcc
CXX_LIB = g++
FC_LIB  = gfortran

# Compile and link flags
CPPFLAGS += -O3 
FFLAGS += -O3 
$FORTRAN_DIR = 
LDFLAGS += $(LIB_DIR)/lib_mlip_cblas.a -L$(FORTRAN_DIR) -lgfortran
LIB_BLAS += $(LIB_DIR)/lib_mlip_cblas.a
CPPFLAGS += -I./cblas

# Extra variables
TARGET_PRERQ = $(LIB_DIR)/lib_mlip_cblas.a
