# the compiler
CC = gcc
CXX = g++

# the compiler flags
CFLAGS = -Wall -g -I. -O3
CXXFLAGS = -Wall -g -I. -O3

LFLAGS = -lm  -g

BIN = bin

MOMENT_OBJ = MomentEq.o
RATE_OBJ = RateEq.o
PAIRS_OBJ = PairEq.o

# the first target. All executables will be made when 'make' is
# invoked with no target
all: moment rate pairs

# general rule how to compile a source file and produce an object file
%.c,%.o:
	$(CXX) $(CFLAGS) -c $<

%.cpp,%.o:
	$(CXX) $(CFLAGS) -c $<

# linking rule for the first executable

moment: $(MOMENT_OBJ)
	$(CXX) $(LFLAGS)  $(MOMENT_OBJ) -o $(BIN)/$@	

rate: $(RATE_OBJ)
	$(CXX) $(LFLAGS)  $(RATE_OBJ) -o $(BIN)/$@

pairs: $(PAIRS_OBJ)
	$(CXX) $(LFLAGS)  $(PAIRS_OBJ) -o $(BIN)/$@

# remove object files and backups
clean:
	rm -f *.o *~ 


PairEq.o: rk_solver.h equation_set.h vm/vec_mat.h
PairEq.o: vm/vm_traits.h vm/vector.h vm/matrix.h vm/vm_types.h vm/vm_math.h
PairEq.o: vm/vm_math.cpp vm/vm_tools.h vm/vm_tools.cpp vm/vm_io.h
PairEq.o: vm/vm_io.cpp ChemicalNetwork.h


MomentEq.o: rk_solver.h equation_set.h vm/vec_mat.h
MomentEq.o: vm/vm_traits.h vm/vector.h vm/matrix.h vm/vm_types.h vm/vm_math.h
MomentEq.o: vm/vm_math.cpp vm/vm_tools.h vm/vm_tools.cpp vm/vm_io.h
MomentEq.o: vm/vm_io.cpp 

RateEq.o: rk_solver.h equation_set.h vm/vec_mat.h 
RateEq.o: vm/vm_traits.h vm/vector.h vm/matrix.h vm/vm_types.h vm/vm_math.h
RateEq.o: vm/vm_math.cpp vm/vm_tools.h vm/vm_tools.cpp vm/vm_io.h
RateEq.o: vm/vm_io.cpp 

