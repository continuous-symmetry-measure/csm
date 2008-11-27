# the compiler
CC = gcc
CXX = g++

# the compiler flags
CFLAGS = -Wall -g -I. 
CXXFLAGS = -Wall -g -I.

LFLAGS = -lm  -g

BIN = bin

MOMENT_OBJ = MomentEq.o
RATE_OBJ = RateEq.o
PAIRS_OBJ = PairEq.o
MASTER_OBJ = MasterEq.o
MASTER_NEWTON_OBJ = MasterNewton.o

# the first target. All executables will be made when 'make' is
# invoked with no target
all: moment rate pairs master

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

master: $(MASTER_OBJ)
	$(CXX) $(LFLAGS)  $(MASTER_OBJ) -o $(BIN)/$@

newton: $(MASTER_NEWTON_OBJ)
	$(CXX) $(LFLAGS)  $(MASTER_NEWTON_OBJ) -o $(BIN)/$@

# remove object files and backups
clean:
	rm -f *.o *~ 
