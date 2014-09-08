/*
 * Author: shadi lahham
 *
 * contains helper functions for the C2 operation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>   // for acos, cos, etc
#include "Molecule.h"
#include "mainhelpers.h"



// ************************************************************
//       global variables
// ************************************************************
int operationOrder = 2;
char* operationName = "INVERSION (S2)";



// ************************************************************
//       function declarations
// ************************************************************

double determinant(double R[3][3]);
void cubic(double a, double b, double c, double *eigen_val);



// ************************************************************
//       function implementations
// ************************************************************

/*
 * creates an antimer from the Molecule. In case of C2 it's just a copy of the
 * original points
 */
double** createAntimer(Molecule* m){

	int i,j;

	 // allocate antimer
	double ** antimer = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		antimer[i] = (double *)malloc(3 * sizeof(double));

	// init antimer - for c2 just copy Molecule pos
	for (i=0;i<m->_size;i++)
		for (j=0;j<3;j++)
			antimer[i][j] = m->_pos[i][j];

	return antimer;
}

/*
 * frees memory allocated for the antimer
 */
void freeAntimer(Molecule* m,double** antimer){

	int i;
    // free antimer
    for (i=0;i<m->_size;i++){
		free(antimer[i]);
	}

	free(antimer);

}

/*
 * s2_estimation in this case of s2
 */
double chiralityFunction(double *parray[], double *anti_parray[],
                        int num_of_atoms, double *dir, double *d_min, double eulerParam[]){

  double csm=MAXDOUBLE;
  int atom;


  *d_min=0.0;
  for(atom = 0; atom < num_of_atoms; ++atom)
  {
    *d_min +=  parray[atom][0]*anti_parray[atom][0];
    *d_min +=  parray[atom][1]*anti_parray[atom][1];
    *d_min +=  parray[atom][2]*anti_parray[atom][2];
  }

  // normalize to 1 and not num of atoms
  //*d_min = (*d_min * -1)/(double)num_of_atoms;
  *d_min = -*d_min;

   csm = (1 - *d_min);
   dir[0] = dir[1] = dir[2] = 0.0;

   return(csm);
}

void foldingUnfolding(double *parray[], double *anti_parray[],double *uarray[],
                      int num_of_atoms, double *dir_cos,double d_min, double eulerParam[]){
  int atom;
  double x, y, z;

  for(atom = 0; atom < num_of_atoms; ++atom)
  {

    x = -(anti_parray[atom][0]);
    y = -(anti_parray[atom][1]);
    z = -(anti_parray[atom][2]);

    uarray[atom][0] = (x+parray[atom][0])/2.0*d_min;
    uarray[atom][1] = (y+parray[atom][1])/2.0*d_min;
    uarray[atom][2] = (z+parray[atom][2])/2.0*d_min;
  }
}
