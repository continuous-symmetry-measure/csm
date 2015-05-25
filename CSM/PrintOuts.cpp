/*
* Various print-out utility functions.
*
* Taken from mainRot.cpp during the 2014 organization.
* Created by Itay Zandbank
*/

#include "PrintOuts.h"

#include <stdio.h>
#include "Molecule.h"
#include <string>
#include "options.h"

/*
* prints the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
/*void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out, double* localCSM)
{

	int i, j;
	printf("%s: %.6lf\n", options.opName.c_str(), fabs(csm));
	fprintf(out, "%s: %.4lf\n", options.opName.c_str(), fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n%i\n", m->size());

	for (i = 0; i<m->size(); i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
			m->symbol(i), m->pos()[i][0], m->pos()[i][1], m->pos()[i][2]);
	}

	for (i = 0; i < m->size(); i++) {
		fprintf(out, "%d ", i + 1);
		for (j = 0; j < m->valency(i); j++) {
			fprintf(out, "%d ", m->adjacent(i, j) + 1);
		}
		fprintf(out, "\n");
	}

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n%i\n", m->size());

	for (i = 0; i<m->size(); i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
			m->symbol(i), outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
	}

	for (i = 0; i < m->size(); i++) {
		fprintf(out, "%d ", i + 1);
		for (j = 0; j < m->valency(i); j++) {
			fprintf(out, "%d ", m->adjacent(i, j) + 1);
		}
		fprintf(out, "\n");
	}

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (options.printNorm) {
		printf("NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf("DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf("NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
	}

	if (options.printLocal) {
		double sum = 0;
		fprintf(out, "\nLocal CSM: \n");
		for (i = 0; i < m->size(); i++) {
			sum += localCSM[i];
			fprintf(out, "%s %7lf\n", m->symbol(i), localCSM[i]);
		}
		fprintf(out, "\nsum: %7lf\n", sum);
	}



}

*/

/*
* prints PDB ATOM tags
*/
void printPDBATOM(Molecule* m, FILE* f, char** sym, double** pos)
{
	int i;
	for (i = 0; i<m->size(); i++){
		fprintf(f, "ATOM  %5d %2s                %8.3lf%8.3lf%8.3lf                      %2s\n",
			i + 1, sym[i], pos[i][0], pos[i][1], pos[i][2], sym[i]);
	}
}

/*
* prints PDB CONECT tags
*/
void printPDBCONNECT(Molecule* m, FILE* f){
	int i, j;
	for (i = 0; i<m->size(); i++){
		fprintf(f, "CONECT%5d", i + 1);
		for (j = 0; j< m->valency(i); j++){
			if ((j>0) && (!(j % 4)))
				fprintf(f, "\nCONECT%5d", i + 1);
			fprintf(f, "%5d", m->adjacent(i, j) + 1);
		}
		fprintf(f, "\n");
	}
}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out)
{

	//	out = stdout;
	// print PDB file
	fprintf(out, "MODEL        1\n");

	printPDBATOM(m, out, m->symbols(), m->pos());

	printPDBCONNECT(m, out);

	fprintf(out, "ENDMDL\n");
	fprintf(out, "MODEL        2\n");

	printPDBATOM(m, out, m->symbols(), outAtoms);

	printPDBCONNECT(m, out);

	fprintf(out, "ENDMDL\n");

	// print results to screen


	if (options.writeOpenu)
		printf("SV* %.4lf *SV\n", fabs(csm));
	else
		printf("%s: %.4lf\n", options.opName.c_str(), fabs(csm));

	if (options.printNorm) {
		printf("NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf("DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf("NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
	}

}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
/*
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, const char *fname, double* localCSM) 
{

	fprintf(out, "%s: %.4lf\n", options.opName.c_str(), fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	// TODO - should we print the centered molecule, or the original one (and, accordingly, the symmetric struct)

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n");

	updateCoordinates(mol, m->pos());

	writeMolecule(mol, options.format, out, fname);

	updateCoordinates(mol, outAtoms);

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n");

	writeMolecule(mol, options.format, out, fname);

	// print results to screen

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (options.writeOpenu)
		printf("SV* %.4lf *SV\n", fabs(csm));
	else
		printf("%s: %.4lf\n", options.opName.c_str(), fabs(csm));

	if (options.printNorm) {
		printf("NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf("DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf("NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
	}

	if (options.printLocal) {
		double sum = 0;
		int i;
		fprintf(out, "\nLocal CSM: \n");
		for (i = 0; i < m->size(); i++) {
			sum += localCSM[i];
			fprintf(out, "%s %7lf\n", m->symbol(i), localCSM[i]);
		}
		fprintf(out, "\nsum: %7lf\n", sum);
	}
}
*/
