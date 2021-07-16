# This is a simple Makefile for the statistic_calc programme compilation on the Okeanos cluster
# If you want to keep all the objects already compiled, remove the "clean_objects" dependency from "all" rule
# The project uses e.x. the std::vector container, therefore a compiler supporting the (at least) c++11 library standard is needed
# Compile with GNU C++ compilers, provided by loading PrgEnv-gnu

# k.topolski2@student.uw.edu.pl

# LOCATIONS OF /lib DIRECTORIES

FFTW_DIR   = /home/topolski/lib
GSL_DIR    = /home/topolski/lib

# LOCATIONS OF /include DIRECTORIES

FFTW_INC       = /home/topolski/include
GSL_INC        = /home/topolski/include

# NECESSARY LIBARIES TO LINK


FFTW_LIB   = -lfftw3_omp -lfftw3 #-lfftw3f
GSL_LIB    = -lgsl -lgslcblas

DIRS    =
INCS    =
LIBS    =

DIRS    += -L $(FFTW_DIR)
DIRS    += -L $(GSL_DIR)

INCS    += -I $(FFTW_INC)
INCS    += -I $(GSL_INC)

LIBS += $(FFTW_LIB) 
LIBS += $(GSL_LIB)

CC	= g++ -std=c++17
CCFLAGS = -lm -fopenmp -lstdc++fs
LDFLAGS = -rdynamic $(INCS) $(DIRS) 


TARGETS := statistics_calc

.PHONY: all clean 

all: $(TARGETS)

clean:
	@rm -f $(TARGETS) 
	@echo "Removing the executable"

$(TARGETS): 
	$(CC) statistics_calc.cpp -o $@ $(LIBS) $^ $(CCFLAGS) $(LDFLAGS)
	


