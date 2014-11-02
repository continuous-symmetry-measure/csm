# the compiler
CC = gcc
CXX = g++

# the compiler flags
CFLAGS = -Wall -O3 -DUSE_ORDER2_PERMUTER -I../../openbabel-files/include/openbabel-2.0 -I../../openbabel-files/os/include -g -I/usr/local/include
CXXFLAGS = -Wall -O3 -DBOOST_ALL_DYN_LINK -DUSE_ORDER2_PERMUTER -I../../openbabel-files/include -I../../openbabel-files/os/include/openbabel-2.0 -I../../include -I/usr/local/include -g -std=c++11
LFLAGS = -L. -L../../openbabel-files/os/lib -L/usr/local/lib -lopenbabel -lboost_log-mt -lboost_log_setup-mt -lboost_thread-mt -lboost_system-mt -lpthread

BIN = .

OBJCSM = Molecule.o parseFunctions.o permuter.o groupPermuter.o elements.o logging.o \
         mainRot.o rpoly.o babelAdapter.o math_wrappers.o \
         drand48.o math_utils.o calculations.o options.o PrintOuts.o

all: csm 

# general rule how to compile a source file and produce an object file
vpath %.cpp ../
vpath %.c ../

%.c,%.o:
	$(CXX) $(CFLAGS) -c $<

%.cpp,%.o:
	$(CXX) $(CFLAGS) -c $<

# linking rule for the first executable
#c2: $(OBJ1) libcsm.a
#	$(CXX) $(LFLAGS) -lcsm $(OBJ1) -o $(BIN)/$@

# linking rule for the second executable
#mir: $(OBJ2) libcsm.a
#	$(CXX) $(LFLAGS) -lcsm $(OBJ2) -o $(BIN)/$@

# linking rule for the third executable
#chiral: $(OBJ3) libcsm.a
#	$(CXX) $(LFLAGS) -lcsm $(OBJ3) -o $(BIN)/$@

# linking rule for the fourth executable
#inv: $(OBJ4) libcsm.a
#	$(CXX) $(LFLAGS) -lcsm $(OBJ4) -o $(BIN)/$@

# linking rule for the csm executable
csm: $(OBJCSM) 
	$(CXX) $(OBJCSM) $(LFLAGS) -o $(BIN)/$@

# linking rule for the csm2d executable
#csm2d: $(OBJCSM2D) libcsm.a
#	$(CXX) $(LFLAGS) $(OBJCSM2D) -o $(BIN)/$@

csm2pdb: $(OBJ5)
	$(CC) $(LFLAGS) $(OBJ5) -o $(BIN)/$@

#libcsm.a: $(LIBOBJ)
#	ar -r $@ $(LIBOBJ)
#	ranlib $@

# make .o file dependencies
depend:
	echo -e "###\n" >> makefile
	$(CXX) -MM *.c >> makefile

# remove object files and backups
clean:
	rm -f *.o *~ libcsm.a

