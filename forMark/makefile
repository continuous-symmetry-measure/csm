# the compiler
CC = gcc

# the compiler flags
CFLAGS = -g -O3 -Wall  -DUSE_ORDER2_PERMUTER 

LFLAGS = -lm -L. 

BIN = bin


OBJCSM =  Molecule.o parseFunctions.o permuter.o groupPermuter.o main.o rpoly.o tqli.o tred2.o nrutil.o c2helpers.o pythag.o icosahedron.c

# the first target. All executables will be made when 'make' is
# invoked with no target
all: cnv

# general rule how to compile a source file and produce an object file
%.c,%.o:
	$(CC) $(CFLAGS) -c $<

# linking rule for the csm executable
cnv: $(OBJCSM)
	$(CC) $(LFLAGS) $(OBJCSM) -o $(BIN)/$@

# make .o file dependencies
depend:
	echo -e "###\n" >> makefile
	$(CC) -MM *.c >> makefile

# remove object files and backups
clean:
	rm -f *.o *~ 

