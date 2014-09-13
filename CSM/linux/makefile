# the compiler
CC = gcc
CXX = g++

# the compiler flags
CFLAGS = -Wall -O3 -DUSE_ORDER2_PERMUTER -I../../openbabel-files/include/openbabel-2.0 -I../../openbabel-files/os/include -g
CXXFLAGS = -Wall -O3 -DUSE_ORDER2_PERMUTER -I../../openbabel-files/include -I../../openbabel-files/os/include/openbabel-2.0 -g -std=c++11
LFLAGS = -L. -L../../openbabel-files/os/lib -lopenbabel

BIN = .

# object files for executables
#LIBOBJ = Molecule.o parseFunctions.o permuter.o groupPermuter.o main.o babelAdapter.o
#OBJ1 = c2helpers.o 
#OBJ2 = mirhelpers.o  
#OBJ3 = chiralhelpers.o rpoly.o nrutil.o
#OBJ4 = invhelpers.o 
#OBJ5 = csm2pdb.o parseFunctions.o Molecule.o

OBJCSM =  Molecule.o parseFunctions.o permuter.o groupPermuter.o mainRot.o rpoly.o tqli.o tred2.o nrutil.o pythag.o babelAdapter.o

#OBJCSM2D = Molecule.o parseFunctions.o permuter.o groupPermuter.o main2d.o rpoly.o tqli.o tred2.o nrutil.o c2helpers.o pythag.o babelAdapter.o




# the first target. All executables will be made when 'make' is
# invoked with no target
#all: libcsm.a c2 mir chiral inv csm csm2pdb csm2d
all: csm 
#csm2pdb

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

