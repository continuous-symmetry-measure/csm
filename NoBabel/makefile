# the compiler
CC = gcc

# the compiler flags
CFLAGS = -g -O3 -Wall  -DUSE_ORDER2_PERMUTER 

LFLAGS = -lm -L. 

BIN = bin

# object files for executables
LIBOBJ = Molecule.o parseFunctions.o permuter.o groupPermuter.o main.o 
OBJ1 = c2helpers.o 
OBJ2 = mirhelpers.o  
OBJ3 = chiralhelpers.o poly_solve.o nrutil.o
OBJ4 = invhelpers.o  

OBJCSM =  Molecule.o parseFunctions.o permuter.o groupPermuter.o mainRot.o rpoly.o tqli.o tred2.o nrutil.o c2helpers.o pythag.o

# the first target. All executables will be made when 'make' is
# invoked with no target
all: libcsm.a c2 mir chiral inv csm

# general rule how to compile a source file and produce an object file
%.c,%.o:
	$(CC) $(CFLAGS) -c $<

# linking rule for the first executable
c2: $(OBJ1)
	$(CC) $(LFLAGS) -lcsm $(OBJ1) -o $(BIN)/$@

# linking rule for the second executable
mir: $(OBJ2)
	$(CC) $(LFLAGS) -lcsm $(OBJ2) -o $(BIN)/$@

# linking rule for the third executable
chiral: $(OBJ3)
	$(CC) $(LFLAGS) -lcsm $(OBJ3) -o $(BIN)/$@

# linking rule for the fourth executable
inv: $(OBJ4)
	$(CC) $(LFLAGS) -lcsm $(OBJ4) -o $(BIN)/$@

# linking rule for the csm executable
csm: $(OBJCSM)
	$(CC) $(LFLAGS) $(OBJCSM) -o $(BIN)/$@

libcsm.a: $(LIBOBJ)
	ar -r $@ $(LIBOBJ)
	ranlib $@

# make .o file dependencies
depend:
	echo -e "###\n" >> makefile
	$(CC) -MM *.c >> makefile

# remove object files and backups
clean:
	rm -f *.o *~ libcsm.a $(BIN)/*

