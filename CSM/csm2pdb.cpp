/*
* Author: shadi lahham, modified by Amir Zayit
*
* Main body that initiates chirality operations.
*
* deals with input output, and the main logic of the high level calculation
*
*/

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> //for strcmp,strlen
#include "Molecule.h"
}

#define TRUE 1
#define FALSE 0

// function declarations

void printOutputPDB(Molecule* m, FILE *out);


// file pointers
FILE* inFile = NULL;
FILE* outFile = NULL;
char *inFileName = NULL;
char *outFileName = NULL;

char opName[100];

void usage(char *op) {
	printf("Usage: %s input_file output_file\n", op);
}

char *getExtension(char *fname) {
	return strrchr(fname,'.') + 1;
}

/*
* parses the command line parameters
*/
void parseInput(int argc, char *argv[]){

	// check number of arguments
	if (argc < 3){
		usage(argv[0]);
		exit(1);
	}

	// try to open infile for reading
	inFileName = argv[1];
	if ((inFile = fopen(inFileName, "rt")) == NULL){		
		printf("Failed to open data file %s\n", inFileName);		
		exit(1);
	}


	// try to open outfile for writing
	outFileName = argv[2];
	if ((outFile = fopen(outFileName, "w")) == NULL){
		printf("Failed to open output file %s for writing\n", outFileName);		
		exit(1);
	}
}

// ************************************************************
//       main program
// ************************************************************

/*
* main funciton - check valid parameters, parse molecule and call chirality Operation
*/
int main(int argc, char *argv[]){

	// init options
	parseInput(argc,argv);

	// try to read molecule from infile
	Molecule* m;
	m = createMolecule(inFile,stdout,FALSE);
	printOutputPDB(m, outFile);
	freeMolecule(m);

	fclose(inFile);
	fclose(outFile);

	return 0;
}
/*
* prints PDB ATOM tags
*/
void printPDBATOM(Molecule* m,FILE* f,char** sym,double** pos){
	int i;
	for(i=0; i<m->_size; i++){
		fprintf(f,"ATOM  %5d %2s                %8.3lf%8.3lf%8.3lf                      %2s\n",
			i+1,sym[i],pos[i][0],pos[i][1],pos[i][2],sym[i]);
	}
}

/*
 * prints PDB CONECT tags
 */
void printPDBCONNECT(Molecule* m,FILE* f){
	int i,j;
	for(i=0; i<m->_size; i++){
		fprintf(f,"CONECT%5d",i +1);
		for ( j=0;  j< m->_valency[i] ; j++ ){
			if ((j>0) && (!(j%4)))
				fprintf(f,"\nCONECT%5d",i +1);
			fprintf(f,"%5d",m->_adjacent[i][j] +1);
		}
		fprintf(f,"\n");
	}
}

/*
 * prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
 */
void printOutputPDB(Molecule* m, FILE *out){

	//	out = stdout;
	// print PDB file
	fprintf(out,"MODEL        1\n");

	printPDBATOM(m,out,m->_symbol,m->_pos);

	printPDBCONNECT(m,out);

	fprintf(out,"ENDMDL\n");
}
