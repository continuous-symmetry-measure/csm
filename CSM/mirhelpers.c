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
char* operationName = "MIRROR SYMMETRY";



// ************************************************************
//       function declarations
// ************************************************************
void calc_ref_plane(int num_of_atoms, double *parray[], double *anti_parray[],
                    double *eigen_val_min, double *eigen_val_max,
                    double dir_cos[]);
double determinant(double R[3][3]);
void cubic(double a, double b, double c, double *eigen_val);



// ************************************************************
//       function implementations
// ************************************************************

/*
 * creates an antimer from the Molecule. In case of Mirror it's just a copy of the
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
 * mirror_estimation in this case of mirror
 */
double chiralityFunction(double *parray[], double *anti_parray[],
                        int num_of_atoms, double *dir, double *d_min, double eulerParam[])

{
  double csm=MAXDOUBLE;
  double eigen_val_min, eigen_val_max;
  int atom;

  calc_ref_plane(num_of_atoms, parray, anti_parray, &eigen_val_min,
		 &eigen_val_max, dir);

  *d_min=0.0;
  for(atom = 0; atom < num_of_atoms; ++atom)
  {
    *d_min +=  parray[atom][0]*anti_parray[atom][0];
    *d_min +=  parray[atom][1]*anti_parray[atom][1];
    *d_min +=  parray[atom][2]*anti_parray[atom][2];
  }

  // Normalize to 1 and not num_of_atoms
  //*d_min = (*d_min-2.0*eigen_val_min)/(double)num_of_atoms;
//  *d_min = (*d_min-2.0*eigen_val_min);

   *d_min = (*d_min-2.0*eigen_val_min);   

   csm = 1 - *d_min;   
   return(csm);
}

void foldingUnfolding(double *parray[], double *anti_parray[],double *uarray[],
                      int num_of_atoms, double *dir_cos,double d_min, double eulerParam[]){

  int atom;
  double mul_scale, x, y, z;

  for(atom = 0; atom < num_of_atoms; ++atom)
  {
    mul_scale = anti_parray[atom][0]*dir_cos[0]+
                anti_parray[atom][1]*dir_cos[1]+
                anti_parray[atom][2]*dir_cos[2];

    x = anti_parray[atom][0]-2.0*mul_scale*dir_cos[0];
    y = anti_parray[atom][1]-2.0*mul_scale*dir_cos[1];
    z = anti_parray[atom][2]-2.0*mul_scale*dir_cos[2];

    uarray[atom][0] = (x+parray[atom][0])/2.0*d_min;
    uarray[atom][1] = (y+parray[atom][1])/2.0*d_min;
    uarray[atom][2] = (z+parray[atom][2])/2.0*d_min;
  }
}

void calc_ref_plane(int num_of_atoms, double *parray[], double *anti_parray[],
                    double *eigen_val_min, double *eigen_val_max,
                    double dir_cos[])
{
  int i, j;
  int atom;
  double a,b,c,tmp1;
  double coef[3][3]={{.0,.0,.0},{.0,.0,.0},{.0,.0,.0}};
  double comp_det[3][3], tmp[3][3];
  double eigen_val[3];


  for(atom = 0; atom < num_of_atoms; ++atom)
     {
       coef[0][0] +=  parray[atom][0]*anti_parray[atom][0];
       coef[1][1] +=  parray[atom][1]*anti_parray[atom][1];
       coef[2][2] +=  parray[atom][2]*anti_parray[atom][2];
       coef[1][0] += (parray[atom][0]*anti_parray[atom][1]+
         parray[atom][1]*anti_parray[atom][0])/2.0;
       coef[2][0] += (parray[atom][0]*anti_parray[atom][2]+
         parray[atom][2]*anti_parray[atom][0])/2.0;
       coef[2][1] += (parray[atom][1]*anti_parray[atom][2]+
         parray[atom][2]*anti_parray[atom][1])/2.0;
     }
     coef[0][1] = coef[1][0];
     coef[0][2] = coef[2][0];
     coef[1][2] = coef[2][1];

     a = -(coef[0][0]+coef[1][1]+coef[2][2]);
     b =  (coef[0][0]*coef[1][1]-coef[1][0]*coef[0][1])+
       (coef[0][0]*coef[2][2]-coef[2][0]*coef[0][2])+
       (coef[1][1]*coef[2][2]-coef[1][2]*coef[2][1]);

     c = -determinant(coef);
			
     cubic(a, b, c, eigen_val);


     *eigen_val_min = 1.0e10;
	 for (i=0; i<3; i++) {
       *eigen_val_min = (*eigen_val_min < eigen_val[i]) ?
       *eigen_val_min : eigen_val[i];
	 }

   for (i=0; i<3; i++)
       for (j=0; j<3; j++)
          tmp[i][j] = coef[i][j];

    for (i=0; i<3; i++)
       tmp[i][i] = coef[i][i] - *eigen_val_min;

    comp_det[0][0] =  (tmp[1][1]*tmp[2][2]-tmp[1][2]*tmp[2][1]);
    comp_det[1][1] =  (tmp[0][0]*tmp[2][2]-tmp[0][2]*tmp[2][0]);
    comp_det[2][2] =  (tmp[0][0]*tmp[1][1]-tmp[0][1]*tmp[1][0]);
    comp_det[0][1] = -(tmp[1][0]*tmp[2][2]-tmp[2][0]*tmp[1][2]);
    comp_det[0][2] =  (tmp[1][0]*tmp[2][1]-tmp[2][0]*tmp[1][1]);
    comp_det[1][2] = -(tmp[0][0]*tmp[2][1]-tmp[2][0]*tmp[0][1]);
    comp_det[1][0] = comp_det[0][1];
    comp_det[2][0] = comp_det[0][2];
    comp_det[2][1] = comp_det[1][2];

    for (j=0; j<3; j++)
       if((fabs(comp_det[0][j]) > 1.0e-10) ||
          (fabs(comp_det[1][j]) > 1.0e-10) ||
          (fabs(comp_det[2][j]) > 1.0e-10))
          {
          dir_cos[0] = comp_det[0][j];
          dir_cos[1] = comp_det[1][j];
          dir_cos[2] = comp_det[2][j];
          break;
          }

    tmp1 = sqrt(dir_cos[0]*dir_cos[0]+
                dir_cos[1]*dir_cos[1]+
                dir_cos[2]*dir_cos[2]);

    for (i=0; i<3; i++)
       dir_cos[i] /= tmp1;
}

double determinant(double R[3][3])
{
  double tmp = 0.0;

  tmp += R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1]);
  tmp += R[0][1]*(R[1][2]*R[2][0]-R[1][0]*R[2][2]);
  tmp += R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
  return(tmp);
}

void cubic(double a, double b, double c, double *eigen_val)
{
	double eps = 1e-10;
	double a2, q, r, z,a3;
	if(fabs(a) < MINDOOUBLE && fabs(b) < MINDOOUBLE && fabs(c) < MINDOOUBLE)
	{
		eigen_val[0] = eigen_val[1] = eigen_val[2] = 0.0;
		return;
	}
	a2 = a*a;
	a3  = a/3.0;
	q = (a2 - 3.0 * b) / 9.0;
	r = (2.0*a*a2 - 9.0*a*b + 27.0*c) / 54.0;
	z = q*q*q-r*r;
	if (fabs(z) < eps) {
		if (fabs(q) < eps) {
			eigen_val[0] = -a3;	
			eigen_val[1] = -a3;
			eigen_val[2] = -a3;
		} else if (r < 0.0) {
			eigen_val[0] = 2.0 * sqrt(q) - a3;
			eigen_val[1] = - sqrt(q) - a3;
			eigen_val[2] = - sqrt(q) - a3;
		} else {
			eigen_val[0] = sqrt(q) - a3;
			eigen_val[1] = -2.0 * sqrt(q) - a3;
			eigen_val[2] = -2.0 * sqrt(q) - a3;
		}
	} else {
		if (z < 0.0) {
			double z2 = pow(sqrt(-z) + fabs(r),1.0/3.0);
			double z1 =	z2 + q/z2;
			if (r > 0) z1 = -z1;
			eigen_val[0] = z1-a3;
			eigen_val[1] = z1-a3;
			eigen_val[2] = z1-a3;
		} else {
			double theta = r/sqrt(q*q*q);
			double q2 = -2.0*sqrt(q);
			if (theta > 0.9999999999) {
				theta = 0.0;
			} else if (theta < -0.9999999999) {
				theta = PI;
			} else {
				theta = acos(theta);
			}
			eigen_val[0] = q2*cos(theta/3.0)-a3;
			eigen_val[1] = q2*cos((theta + 2.0 * PI)/3.0)-a3;
			eigen_val[2] = q2*cos((theta + 4.0 * PI)/3.0)-a3;
		}
	}
}

