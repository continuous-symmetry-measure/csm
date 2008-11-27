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
char* operationName = "CHIRALITY MEASURE";



// ************************************************************
//       function declarations
// ************************************************************
void calc_ref_plane(int num_of_atoms, double *parray[], double *anti_parray[],
                    double *eigen_val_min, double *eigen_val_max,
                    double dir_cos[], double euler_param[]);
double determinant(double R[3][3]);

// in poly_solve
void zrhqr(double a[], int m, double rtr[], double rti[]);

// ************************************************************
//       function implementations
// ************************************************************

/*
 * creates an antimer from the Molecule. In case of chrality the first
 * coordinant x is multiplied by (-1)
 */
double** createAntimer(Molecule* m){

	int i;

	 // allocate antimer
	double ** antimer = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		antimer[i] = (double *)malloc(3 * sizeof(double));

	// init antimer - for c2 just copy Molecule pos
	for (i=0;i<m->_size;i++){
		antimer[i][0] = - m->_pos[i][0];
		antimer[i][1] =   m->_pos[i][1];
		antimer[i][2] =   m->_pos[i][2];
	}

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
 * (originally called) cn_estimation in this case of chirality.
 */
double chiralityFunction(double *parray[], double *anti_parray[],
                        int num_of_atoms, double *dir, double *d_min, double eulerParam[])
{
  double csm=MAXDOUBLE;
  double eigen_val_min, eigen_val_max;

  calc_ref_plane(num_of_atoms, parray, anti_parray, &eigen_val_min,
		 &eigen_val_max, dir, eulerParam);

  // normalize to 1 and not num_of_atoms
  //*d_min = eigen_val_max/(num_of_atoms);
  *d_min = eigen_val_max;
  csm = 1 - *d_min;

   return(csm);
}

/*
 * Normalizes uarray - helper for dina correction in folding unfolding
 * not an elegant solution! Should make outAtoms/uarray a molecule instance
 * returns one [TRUE] if successful, zero[FALSE] otherwise
 */
int normalizePositions(double *uarray[], int size){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	for(i=0; i< size; i++){
		x_avg += uarray[i][0];
		y_avg += uarray[i][1];
		z_avg += uarray[i][2];
	}
	x_avg /= (double)(size);
	y_avg /= (double)(size);
	z_avg /= (double)(size);

	norm = 0.0;
	for(i=0; i< size; i++){
		tmp = SQR(uarray[i][0]-x_avg) +
		      SQR(uarray[i][1]-y_avg) +
		      SQR(uarray[i][2]-z_avg);
		norm += tmp;
	}
	//normalize to 1 and not num of atoms
	// norm = sqrt(norm/(double)size);
	norm = sqrt(norm);

	if(norm < MINDOOUBLE)
		return FALSE;

	for(i=0; i< size; i++){
		uarray[i][0] = ((uarray[i][0] - x_avg) / norm);
		uarray[i][1] = ((uarray[i][1] - y_avg) / norm);
		uarray[i][2] = ((uarray[i][2] - z_avg) / norm);
	}

	return TRUE;
}

void foldingUnfolding(double *parray[], double *anti_parray[],double *uarray[],
                      int num_of_atoms, double *dir_cos,double d_min, double eulerParam[])
{
  int k;
  double x, y, z;
  double rot_matrx[3][3];
  double l,l2,m,m2,p,p2,v,v2;

  l=eulerParam[0]; m=eulerParam[1]; v=eulerParam[2] ;p=eulerParam[3];
  l2=l*l;m2=m*m;v2=v*v;p2=p*p;

  rot_matrx[0][0]=l2-m2-v2+p2; rot_matrx[0][1]=2*(l*m-v*p); rot_matrx[0][2]=2*(v*l+m*p);
  rot_matrx[1][0]=2*(l*m+v*p); rot_matrx[1][1]=m2-v2-l2+p2; rot_matrx[1][2]=2*(m*v-l*p);
  rot_matrx[2][0]=2*(v*l-m*p); rot_matrx[2][1]=2*(m*v+l*p); rot_matrx[2][2]=v2-l2-m2+p2;

  for(k=0; k < num_of_atoms; ++k)
    uarray[k] = (double *)calloc(3,sizeof(double));

  for(k=0; k < num_of_atoms; ++k)
    {
    // ssl - changed from parray to anti_parray
      x=anti_parray[k][0];
      y=anti_parray[k][1];
      z=anti_parray[k][2];

      uarray[k][0] =  rot_matrx[0][0]*x+rot_matrx[1][0]*y+rot_matrx[2][0]*z;
      uarray[k][1] =  rot_matrx[0][1]*x+rot_matrx[1][1]*y+rot_matrx[2][1]*z;
      uarray[k][2] =  rot_matrx[0][2]*x+rot_matrx[1][2]*y+rot_matrx[2][2]*z;

      uarray[k][0] *= d_min;
      uarray[k][1] *= d_min;
      uarray[k][2] *= d_min;
  }

    // ssl - corrected accoding to dina (result = average of original and 'normalized' result)
	// ps: foldingUnfolding is not the best place to do this, try to move it somewhere more logical?

	normalizePositions(uarray,num_of_atoms); //ignore TRUE/FALSE return value - ugly

	for(k=0; k < num_of_atoms; ++k) {
		// calc average
		uarray[k][0] = (uarray[k][0] + parray[k][0]) / 2.0;
		uarray[k][1] = (uarray[k][1] + parray[k][1]) / 2.0;
		uarray[k][2] = (uarray[k][2] + parray[k][2]) / 2.0;
    }
}

void calc_ref_plane(int num_of_atoms, double *parray[], double *anti_parray[],
                    double *eigen_val_min, double *eigen_val_max,
                    double dir_cos[], double euler_param[])
{
  int i, j, k, k_max = -1, k_min = -1;
  int atom;
  double sum2, sin_t2;
  double coef[4][4]={{.0,.0,.0,.0},{.0,.0,.0,.0},{.0,.0,.0,.0},{.0,.0,.0,.0}};
  double comp_det[4][4], tmp[4][4], tmp2[3][3];
  double eigen_val[4];
  double rtr[5], rti[5];
  double poly_c[5];


  for(atom = 0; atom < num_of_atoms; ++atom)
     {
      coef[0][0] +=  ( parray[atom][0]*anti_parray[atom][0]-
                       parray[atom][1]*anti_parray[atom][1]-
                       parray[atom][2]*anti_parray[atom][2]);
      coef[1][1] +=  (-parray[atom][0]*anti_parray[atom][0]+
                       parray[atom][1]*anti_parray[atom][1]-
                       parray[atom][2]*anti_parray[atom][2]);
      coef[2][2] +=  (-parray[atom][0]*anti_parray[atom][0]-
                       parray[atom][1]*anti_parray[atom][1]+
                       parray[atom][2]*anti_parray[atom][2]);
      coef[3][3] +=  ( parray[atom][0]*anti_parray[atom][0]+
                       parray[atom][1]*anti_parray[atom][1]+
                       parray[atom][2]*anti_parray[atom][2]);

      coef[1][0] += (parray[atom][0]*anti_parray[atom][1]+
                     parray[atom][1]*anti_parray[atom][0]);
      coef[2][0] += (parray[atom][0]*anti_parray[atom][2]+
                     parray[atom][2]*anti_parray[atom][0]);
      coef[3][0] += (parray[atom][1]*anti_parray[atom][2]-
                     parray[atom][2]*anti_parray[atom][1]);
      coef[2][1] += (parray[atom][1]*anti_parray[atom][2]+
                     parray[atom][2]*anti_parray[atom][1]);
      coef[3][1] += (parray[atom][2]*anti_parray[atom][0]-
                     parray[atom][0]*anti_parray[atom][2]);
      coef[3][2] += (parray[atom][0]*anti_parray[atom][1]-
                     parray[atom][1]*anti_parray[atom][0]);

     }
  coef[0][1] = coef[1][0];
  coef[0][2] = coef[2][0];
  coef[0][3] = coef[3][0];
  coef[1][2] = coef[2][1];
  coef[1][3] = coef[3][1];
  coef[2][3] = coef[3][2];

  poly_c[4] = 1.0;
  poly_c[3] = 0.0;
  poly_c[2] = -coef[0][3]*coef[3][0]-coef[0][1]*coef[1][0]+
               coef[0][0]*coef[1][1]+coef[0][0]*coef[2][2]-
               coef[0][2]*coef[2][0]+coef[0][0]*coef[3][3]-
               coef[1][3]*coef[3][1]-coef[1][2]*coef[2][1]+
               coef[1][1]*coef[2][2]+coef[1][1]*coef[3][3]-
               coef[2][3]*coef[3][2]+coef[2][2]*coef[3][3];
  poly_c[1] = -coef[0][2]*coef[3][0]*coef[2][3]+
               coef[0][2]*coef[1][1]*coef[2][0]-
               coef[0][0]*coef[1][1]*coef[2][2]-
               coef[0][2]*coef[1][0]*coef[2][1]+
               coef[0][2]*coef[2][0]*coef[3][3]-
               coef[1][3]*coef[2][1]*coef[3][2]+
               coef[0][1]*coef[1][0]*coef[3][3]+
               coef[1][2]*coef[2][1]*coef[3][3]-
               coef[0][3]*coef[1][0]*coef[3][1]-
               coef[0][0]*coef[1][1]*coef[3][3]-
               coef[0][1]*coef[1][2]*coef[2][0]+
               coef[0][1]*coef[1][0]*coef[2][2]+
               coef[0][0]*coef[1][3]*coef[3][1]+
               coef[0][0]*coef[1][2]*coef[2][1]-
               coef[0][1]*coef[1][3]*coef[3][0]+
               coef[1][1]*coef[3][2]*coef[2][3]-
               coef[0][3]*coef[2][0]*coef[3][2]+
               coef[0][3]*coef[3][0]*coef[2][2]-
               coef[1][2]*coef[3][1]*coef[2][3]+
               coef[0][3]*coef[1][1]*coef[3][0]+
               coef[1][3]*coef[3][1]*coef[2][2]-
               coef[0][0]*coef[2][2]*coef[3][3]+
               coef[0][0]*coef[3][2]*coef[2][3]-
               coef[1][1]*coef[2][2]*coef[3][3];

  poly_c[0] =  coef[0][0]*coef[1][1]*coef[2][2]*coef[3][3]-
               coef[0][0]*coef[1][1]*coef[2][3]*coef[3][2]-
               coef[0][0]*coef[1][2]*coef[2][1]*coef[3][3]+
               coef[0][0]*coef[1][2]*coef[3][1]*coef[2][3]+
               coef[0][0]*coef[1][3]*coef[2][1]*coef[3][2]-
               coef[0][0]*coef[1][3]*coef[3][1]*coef[2][2]-
               coef[0][1]*coef[1][0]*coef[2][2]*coef[3][3]+
               coef[0][1]*coef[1][0]*coef[3][2]*coef[2][3]+
               coef[0][1]*coef[1][2]*coef[2][0]*coef[3][3]-
               coef[0][1]*coef[1][2]*coef[3][0]*coef[2][3]-
               coef[0][1]*coef[1][3]*coef[2][0]*coef[3][2]+
               coef[0][1]*coef[1][3]*coef[3][0]*coef[2][2]+
               coef[0][2]*coef[1][0]*coef[2][1]*coef[3][3]-
               coef[0][2]*coef[1][0]*coef[3][1]*coef[2][3]-
               coef[0][2]*coef[1][1]*coef[2][0]*coef[3][3]+
               coef[0][2]*coef[1][1]*coef[3][0]*coef[2][3]+
               coef[0][2]*coef[1][3]*coef[2][0]*coef[3][1]-
               coef[0][2]*coef[1][3]*coef[3][0]*coef[2][1]-
               coef[0][3]*coef[1][0]*coef[2][1]*coef[3][2]+
               coef[0][3]*coef[1][0]*coef[3][1]*coef[2][2]+
               coef[0][3]*coef[1][1]*coef[2][0]*coef[3][2]-
               coef[0][3]*coef[1][1]*coef[3][0]*coef[2][2]-
               coef[0][3]*coef[1][2]*coef[2][0]*coef[3][1]+
               coef[0][3]*coef[1][2]*coef[3][0]*coef[2][1];

  zrhqr(poly_c, 4, rtr, rti);

  for(k=0; k<4; ++k)
     eigen_val[k] = rtr[k+1];

  *eigen_val_max = -MAXDOUBLE;
  *eigen_val_min =  MAXDOUBLE;
  for(k=0; k<4; ++k)
     {
     if(eigen_val[k] > *eigen_val_max)
        { *eigen_val_max = eigen_val[k]; k_max = k; }
     if(eigen_val[k] < *eigen_val_min)
        { *eigen_val_min = eigen_val[k]; k_min = k; }
     }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      tmp[i][j] = coef[i][j];

  for(i=0; i<4; i++)
     tmp[i][i] = coef[i][i] - eigen_val[k_max];

  tmp2[0][0] = tmp[1][1]; tmp2[0][1] = tmp[1][2]; tmp2[0][2] = tmp[1][3];
  tmp2[1][0] = tmp[2][1]; tmp2[1][1] = tmp[2][2]; tmp2[1][2] = tmp[2][3];
  tmp2[2][0] = tmp[3][1]; tmp2[2][1] = tmp[3][2]; tmp2[2][2] = tmp[3][3];
  comp_det[0][0] = determinant(tmp2);

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][2]; tmp2[0][2] = tmp[0][3];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][2]; tmp2[1][2] = tmp[2][3];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][2]; tmp2[2][2] = tmp[3][3];
  comp_det[1][1] = determinant(tmp2);

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][1]; tmp2[0][2] = tmp[0][3];
  tmp2[1][0] = tmp[1][0]; tmp2[1][1] = tmp[1][1]; tmp2[1][2] = tmp[1][3];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][3];
  comp_det[2][2] = determinant(tmp2);

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][1]; tmp2[0][2] = tmp[0][2];
  tmp2[1][0] = tmp[1][0]; tmp2[1][1] = tmp[1][1]; tmp2[1][2] = tmp[1][2];
  tmp2[2][0] = tmp[2][0]; tmp2[2][1] = tmp[2][1]; tmp2[2][2] = tmp[2][2];
  comp_det[3][3] = determinant(tmp2);

  tmp2[0][0] = tmp[1][0]; tmp2[0][1] = tmp[1][2]; tmp2[0][2] = tmp[1][3];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][2]; tmp2[1][2] = tmp[2][3];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][2]; tmp2[2][2] = tmp[3][3];
  comp_det[0][1] = -determinant(tmp2);
  comp_det[1][0] = comp_det[0][1];

  tmp2[0][0] = tmp[1][0]; tmp2[0][1] = tmp[1][1]; tmp2[0][2] = tmp[1][3];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][1]; tmp2[1][2] = tmp[2][3];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][3];
  comp_det[0][2] = determinant(tmp2);
  comp_det[2][0] = comp_det[0][2];

  tmp2[0][0] = tmp[1][0]; tmp2[0][1] = tmp[1][1]; tmp2[0][2] = tmp[1][2];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][1]; tmp2[1][2] = tmp[2][2];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][2];
  comp_det[0][3] = -determinant(tmp2);
  comp_det[3][0] = comp_det[0][3];

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][1]; tmp2[0][2] = tmp[0][3];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][1]; tmp2[1][2] = tmp[2][3];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][3];
  comp_det[1][2] = -determinant(tmp2);
  comp_det[2][1] = comp_det[1][2];

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][1]; tmp2[0][2] = tmp[0][2];
  tmp2[1][0] = tmp[2][0]; tmp2[1][1] = tmp[2][1]; tmp2[1][2] = tmp[2][2];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][2];
  comp_det[1][3] = determinant(tmp2);
  comp_det[3][1] = comp_det[1][3];

  tmp2[0][0] = tmp[0][0]; tmp2[0][1] = tmp[0][1]; tmp2[0][2] = tmp[0][2];
  tmp2[1][0] = tmp[1][0]; tmp2[1][1] = tmp[1][1]; tmp2[1][2] = tmp[1][2];
  tmp2[2][0] = tmp[3][0]; tmp2[2][1] = tmp[3][1]; tmp2[2][2] = tmp[3][2];
  comp_det[2][3] = -determinant(tmp2);
  comp_det[3][2] = comp_det[2][3];


  for (k=0; k<4; k++)
    if((fabs(comp_det[0][k]) > 1.0e-10) || (fabs(comp_det[1][k]) > 1.0e-10)||
       (fabs(comp_det[2][k]) > 1.0e-10) || (fabs(comp_det[3][k]) > 1.0e-10))
       {
       sum2 = 0.0;
       for (j=0; j<4; j++)
          sum2 += comp_det[j][k]*comp_det[j][k];

       for (j=0; j<4; j++)
          euler_param[j] = comp_det[j][k]/sqrt(sum2);
       break;
       }

  sin_t2 = sqrt(1.0-euler_param[3]*euler_param[3]);
  for (j=0; j<3; j++)
     dir_cos[j] = euler_param[j]/sin_t2;

}

double determinant(double R[][3])
{
  double tmp = 0.0;

  tmp += R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1]);
  tmp += R[0][1]*(R[1][2]*R[2][0]-R[1][0]*R[2][2]);
  tmp += R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
  return(tmp);
}
