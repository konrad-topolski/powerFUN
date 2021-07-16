# This is a simple Makefile for the powerFUN programme compilation on the Okeanos cluster

# Compile in the bash shell, with the library path variable knowing of the shared objects location (can be set in .bashrc)

# If you want to keep all the objects already compiled, remove the "clean_objects" dependency from "all" rule
# The project uses e.x. the std::vector container, therefore a compiler supporting the (at least) c++11 library standard is needed
# WARNING: compile with GNU C++ compilers. There is a bug with multihreaded calculation of density average while using the Intel compiler (icpc) and 
# I have not looked into it yet 

# k.topolski2@student.uw.edu.pl

CC	:= g++ -std=c++17
CCFLAGS := -I /lustre/tetyda/home/topolski/include -L /lustre/tetyda/home/topolski/lib -lfftw3_omp -lfftw3 -lboost_program_options -lm -fopenmp -lstdc++fs
LDFLAGS := -rdynamic 
LIBS    :=


TARGETS := main
MAINS	:= $(addsuffix .o, $(TARGETS) )
SOURCES := $(wildcard *.cpp)
OBJ	:= $(patsubst %.cpp, %.o, $(SOURCES))
DEPS	:= $(wildcard *.h)

.PHONY: all clean rename clean_objects

all: $(TARGETS) rename clean_objects

clean:
	@rm -f $(TARGETS) $(OBJ)
	@echo "Removing all the files Make constructed"
$(OBJ): %.o : %.cpp $(DEPS)
	@echo "Compiling the source files, obtaining object files"
	@$(CC) -c -o $@ $< $(CCFLAGS) $(LDFLAGS) $(LIBS)


$(TARGETS): % : $(filter-out $(MAINS), $(OBJ)) %.o 
	@echo "Linking the objects to make the executable"
	@$(CC) -o $@ $(LIBS) $^ $(CCFLAGS) $(LDFLAGS)

	
clean_objects: $(TARGETS)
	@rm -f $(OBJ)
	@echo "Removing all the object files"
rename: $(TARGETS)
	@ mv main powerFUN
	@echo "Renaming the executable to powerFUN"













