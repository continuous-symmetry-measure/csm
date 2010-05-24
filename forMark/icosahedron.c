/* icosahedron.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static integer c__5 = 5;
static integer c__10 = 10;
static integer c__15 = 15;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!       THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE        !!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*    Written by Max Tegmark, Max-Planck-Institut fuer Physik, Munich */
/*    April 1996 */
/*    Currently I'm at Univ. of Pennsylvania, max@physics.upenn.edu */

/*  WHAT IS IT? */
/*    This FORTRAN package lets the user pixelize the sphere at */
/*    a wide range of resolutions. It was written primarily for */
/*    map-making in astronomy and cosmology. It is also useful */
/*    for doing integrals over the sphere when the integrand is */
/*    expensive to evaluate, so that one wishes to minimize the */
/*    number of points used. */

/*  DOCUMENTATION: */
/*    The package and its purpose is described in detail in */
/*    a postscript file available from */
/*      http://www.mpa-garching.mpg.de/~max/icosahedron.html */
/*    (faster from Europe) and from from */
/*      http://sns.ias.edu.edu/~max/icosahedron.html */
/*    (faster from the US). This site also contains the latest */
/*    version of the source code. */

/*  RULES: */
/*    The package is public domain, which means that you are */
/*    allowed to use it for any non-commercial purpose whatsoever */
/*    free of charge. The only requirement is that you include an */
/*    appropriate acknowledgement in any publications based on */
/*    work where the package has been used. Also, if you */
/*    redistribute the package, you must not remove this text. */

/*  HOW IT WORKS: */
/*    As a supplement to the above-mentioned postscript file, */
/*    here is a brief summary of the nitty-gritty details. */
/*    To use the package, first call the subroutines */
/*    compute_matrices and compute_corners once and for all, */
/*    as in the demo routine below. This precomputes some */
/*    geometric stuff to save time later. */
/*    The RESOLUTION is a positive integer 1, 2, 3, ... that you */
/*    can choose freely. It determines the number of pixels, which */
/*    is given by */
/*       N = 40*resolution*(resolution-1)+12. */
/*    For instance, resolution=13 gives N=6252 pixels. */
/*    The only subroutines you ever need to use are these two: */
/*    * vector2pixel takes a point on the sphere (specified by a */
/*      unit vector) and returns the number of the pixel to which */
/*      this point belongs, i.e., an integer between 0 and N-1. */
/*    * pixel2vector does the opposite: it takes a pixel number and */
/*      computes the corresponding unit vector. */
/*    The subroutine "demo" below illustrates how the routines are used. */
/*    It produces a text file with the coordinates of all pixels */
/*    together with their pixel numbers. It also outputs the pixel */
/*    number reconstructed with vector2pixel, so you can verify that */
/*    it all works as it should by checking that the first two columns */
/*    in the file are identical. */

/*  YOU DON'T NEED TO KNOW THIS: */
/*    The resolution is defined so that the number of pixels along */
/*    the side of a triangular face is 2*resolution, so there are */
/*    resolution*(2*resolution+1) pixels on each of the 20 faces of */
/*    the icosahedron. To avoid counting pixels on edges more than */
/*    once, the edge pixels are split up half-and-half between the */
/*    two faces to which the edge belongs, so each face in fact */
/*    only contains 2*resolution*(resolution-1) pixels if you ignore the corners. */
/*    The 12 corner pixels aren't lumped in with any face at all, */
/*    so you can see them listed separately as the last 12 pixels */
/*    in test.dat if you run the demo. */
/*    This makes 40*resolution*(resolution-1) + 12 pixels all in all. */
/*    Thanks to Christopher Weth for catching typos in an earlier version */
/*    of this documentation! */

/*  FEEDBACK: */
/*    If you have any questions, comments or suggestions regarding */
/*    this package, please email them to me at max@ias.edu. */
/*    Thanks, */
/*    ;-) */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ****************************************************************** */
/* Subroutine */ int icosahedron_(integer *resolution, integer *n, doublereal 
	*x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int compute_corners__(doublereal *);
    static integer i__;
    extern /* Subroutine */ int compute_matrices__(doublereal *);
    static doublereal r__[180]	/* was [20][3][3] */, v[36]	/* was [12][3]
	     */;
    extern /* Subroutine */ int pixel2vector_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), vector2pixel_(
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer pixel;
    static doublereal vector[3];

/* 	resolution = 1	!    12 pixels in total. */
/* 	resolution = 2	!    92 pixels in total. */
/* 	resolution = 3	!   252 pixels in total. */
/* 	resolution = 4	!   492 pixels in total. */
/* 	resolution = 5	!   812 pixels in total. */
/* 	resolution = 6	!  1212 pixels in total. */
/* 	resolution = 7	!  1692 pixels in total. */
/* 	resolution = 8	!  2252 pixels in total. */
/* 	resolution = 9	!  2892 pixels in total. */
/* 	resolution = 10	!  3612 pixels in total. */
/* 	resolution = 11	!  4412 pixels in total. */
/* 	resolution = 12	!  5292 pixels in total. */
/* 	resolution = 13	!  6252 pixels in total. */
/* 	resolution = 14	!  7292 pixels in total. */
/* 	resolution = 15	!  8412 pixels in total. */
/* 	resolution = 16	!  9612 pixels in total. */
    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    *n = (*resolution << 1) * (*resolution - 1);
    *n = *n * 20 + 12;
    compute_matrices__(r__);
    compute_corners__(v);
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	pixel2vector_(&i__, resolution, r__, v, vector);
	vector2pixel_(vector, resolution, r__, v, &pixel);
	x[i__ + 1] = vector[0];
	y[i__ + 1] = vector[1];
	z__[i__ + 1] = vector[2];
    }
    return 0;
} /* icosahedron_ */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!     THESE SUBROUTINES ARE ALL YOU NEED TO CALL FROM OUTSIDE    !!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!	These subroutines convert between unit vectors and         !!! */
/* !!     pixel numbers, and are the only ones that the user of this !!! */
/* !!     package calls repeatedly:				   !!! */
/* !!	  subroutine vector2pixel(vector,resolution,R,v,pixel)     !!! */
/* !!	  subroutine pixel2vector(pixel,resolution,R,v,vector)     !!! */
/* !!                                                                !!! */
/* !!	These subroutines are called only once, in the beginning,  !!! */
/* !!     and compute the necessary rotation matrices and corner     !!! */
/* !!     vectors once and for all:                                  !!! */
/* !!	  subroutine compute_matrices(R)                           !!! */
/* !!	  subroutine compute_corners(v)                            !!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ****************************************************************** */
/* Subroutine */ int vector2pixel_(doublereal *vector, integer *resolution, 
	doublereal *r__, doublereal *v, integer *pixel)
{
    /* Builtin functions */
    /* Subroutine */ /*int s_paus(char *, ftnlen);*/

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    extern /* Subroutine */ int find_face__(doublereal *, doublereal *, 
	    integer *);
    static doublereal x, y;
    extern /* Subroutine */ int find_another_face__(doublereal *, doublereal *
	    , integer *), getmatrix_(integer *, doublereal *, doublereal *), 
	    vecmatmul2_(doublereal *, doublereal *, doublereal *), 
	    tangentplanepixel_(integer *, doublereal *, doublereal *, integer 
	    *, integer *);
    static integer pixperface;
    static doublereal vec[3];
    static integer pix;
    extern /* Subroutine */ int find_corner__(doublereal *, doublereal *, 
	    integer *);
    static integer face, ifail;
    extern /* Subroutine */ int adjust_(doublereal *, doublereal *);

    /* Parameter adjustments */
    v -= 12;
    r__ -= 80;
    --vector;

    /* Function Body */
    
    
    /*
    if (*resolution < 1) {
	s_paus("Resolution must exceed 0", (ftnlen)24);
    }
    */
    
    
    pixperface = (*resolution << 1) * (*resolution - 1);
    find_face__(&vector[1], &r__[80], &face);
    getmatrix_(&face, &r__[80], a);
    vecmatmul2_(a, &vector[1], vec);
    x = vec[0] / vec[2];
    y = vec[1] / vec[2];
    adjust_(&x, &y);
    tangentplanepixel_(resolution, &x, &y, &pix, &ifail);
    if (ifail > 0) {
/* Try the runner-up face: */
	find_another_face__(&vector[1], &r__[80], &face);
	getmatrix_(&face, &r__[80], a);
	vecmatmul2_(a, &vector[1], vec);
	x = vec[0] / vec[2];
	y = vec[1] / vec[2];
	adjust_(&x, &y);
	tangentplanepixel_(resolution, &x, &y, &pix, &ifail);
    }
    *pixel = face * pixperface + pix;
    if (ifail > 0) {
/* The pixel wasn't in any of those two faces, */
/* so it must be a corner pixel. */
	find_corner__(&vector[1], &v[12], &pix);
	*pixel = pixperface * 20 + pix;
    }
    return 0;
} /* vector2pixel_ */

/* ****************************************************************** */
/* Subroutine */ int pixel2vector_(integer *pixel, integer *resolution, 
	doublereal *r__, doublereal *v, doublereal *vector)
{
    /* Builtin functions */
    /* Subroutine */ /*int s_paus(char *, ftnlen);*/
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int unadjust_(doublereal *, doublereal *);
    static doublereal a[9]	/* was [3][3] */, x, y;
    extern /* Subroutine */ int getmatrix_(integer *, doublereal *, 
	    doublereal *), vecmatmul1_(doublereal *, doublereal *, doublereal 
	    *);
    static integer pixperface;
    extern /* Subroutine */ int tangentplanevector_(integer *, integer *, 
	    doublereal *, doublereal *);
    static integer pix, face;
    static doublereal norm;

/* Returns a unit vector pointing towards pixel. */
/* Resolution must be an even, positive integer. */
    /* Parameter adjustments */
    --vector;
    v -= 12;
    r__ -= 80;

    /* Function Body */
    
    /*
    if (*resolution < 1) {
	s_paus("Resolution must exceed 0", (ftnlen)24);
    }
    */
    
    
    pixperface = (*resolution << 1) * (*resolution - 1);
    
    
    /*
    if (*pixel < 0) {
	s_paus("Error: negative pixel number", (ftnlen)28);
    }
    */
    
    /*
    if (*pixel >= pixperface * 20 + 12) {
	s_paus("Error: pixel number too large", (ftnlen)29);
    }
    */
    
    
    if (pixperface > 0) {
	face = *pixel / pixperface;
	if (face > 20) {
	    face = 20;
	}
    } else {
/* There are no pixels at all on the faces - it's all just corners. */
	face = 20;
    }
    pix = *pixel - face * pixperface;
    if (face < 20) {
/* The pixel is on one of the 20 faces: */
	tangentplanevector_(&pix, resolution, &x, &y);
	unadjust_(&x, &y);
	norm = sqrt(x * x + y * y + 1.);
	vector[1] = x / norm;
	vector[2] = y / norm;
	vector[3] = 1. / norm;
	getmatrix_(&face, &r__[80], a);
	vecmatmul1_(a, &vector[1], &vector[1]);
    } else {
/* This is a corner pixel: */


/*
	if (pix > 11) {
	    s_paus("Error: pixel number too big", (ftnlen)27);
	}
	
*/	
	
	vector[1] = v[pix + 12];
	vector[2] = v[pix + 24];
	vector[3] = v[pix + 36];
    }
    return 0;
} /* pixel2vector_ */

/* ****************************************************************** */
/* Subroutine */ int compute_matrices__(doublereal *r__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */, b[9]	/* was [3][3] */, c__[
	    9]	/* was [3][3] */, d__[9]	/* was [3][3] */, e[9]	/* 
	    was [3][3] */;
    static integer i__, j, n;
    static doublereal x;
    extern /* Subroutine */ int getmatrix_(integer *, doublereal *, 
	    doublereal *), putmatrix_(integer *, doublereal *, doublereal *);
    static doublereal cs, ct, pi, sn;
    extern /* Subroutine */ int matmul1_(doublereal *, doublereal *, 
	    doublereal *), matmul2_(doublereal *, doublereal *, doublereal *);

/* On exit, R will contain the 20 rotation matrices */
/* that rotate the 20 icosahedron faces */
/* into the tangent plane */
/* (the horizontal plane with z=1). */
/* Only called once, so speed is irrelevant. */
    /* Parameter adjustments */
    r__ -= 80;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    a[i__ + j * 3 - 4] = 0.;
	    b[i__ + j * 3 - 4] = 0.;
	    c__[i__ + j * 3 - 4] = 0.;
	    d__[i__ + j * 3 - 4] = 0.;
	}
    }
    pi = 3.14159265358979;
    x = pi * 2. / 5.;
    cs = cos(x);
    sn = sin(x);
    a[0] = cs;
    a[3] = -sn;
    a[1] = sn;
    a[4] = cs;
    a[8] = 1.;
/* A rotates by 72 degrees around the z-axis. */
    x = pi / 5.;
    ct = cos(x) / sin(x);
    cs = ct / sqrt(3.);
    sn = sqrt(1 - ct * ct / 3.);
    c__[0] = 1.;
    c__[4] = cs;
    c__[7] = -sn;
    c__[5] = sn;
    c__[8] = cs;
/* C rotates around the x-axis so that the north pole */
/* ends up at the center of face 1. */
    cs = -.5f;
    sn = sqrt(3.) / 2.;
    d__[0] = cs;
    d__[3] = -sn;
    d__[1] = sn;
    d__[4] = cs;
    d__[8] = 1.;
/* D rotates by 120 degrees around z-axis. */
    matmul1_(c__, d__, e);
    matmul2_(e, c__, b);
/* B rotates face 1 by 120 degrees. */
/* B = CDC^t */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    e[i__ + j * 3 - 4] = 0.;
	}
	e[i__ + i__ * 3 - 4] = 1.;
    }
/* Now E is the identity matrix. */
    putmatrix_(&c__0, &r__[80], e);
    matmul1_(b, a, e);
    matmul1_(b, e, e);
    putmatrix_(&c__5, &r__[80], e);
    matmul1_(e, a, e);
    putmatrix_(&c__10, &r__[80], e);
    matmul1_(e, b, e);
    matmul1_(e, b, e);
    matmul1_(e, a, e);
    putmatrix_(&c__15, &r__[80], e);
    for (n = 0; n <= 15; n += 5) {
	getmatrix_(&n, &r__[80], e);
	for (i__ = 1; i__ <= 4; ++i__) {
	    matmul1_(a, e, e);
	    i__1 = n + i__;
	    putmatrix_(&i__1, &r__[80], e);
	}
    }
/* Now the nth matrix in R will rotate */
/* face 1 into face n. */
/* Multiply by C so that they will rotate */
/* the tangent plane into face n instead: */
    for (n = 0; n <= 19; ++n) {
	getmatrix_(&n, &r__[80], e);
	matmul1_(e, c__, e);
	putmatrix_(&n, &r__[80], e);
    }
    return 0;
} /* compute_matrices__ */

/* ****************************************************************** */
/* Subroutine */ int compute_corners__(doublereal *v)
{
    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), sin(doublereal), cos(
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal z__, pi, rho, dphi;

/* On exit, v will contain unit vectors pointing toward */
/* the 12 icoshedron corners. */
    /* Parameter adjustments */
    v -= 12;

    /* Function Body */
    pi = atan(1.) * 4.;
    dphi = pi * 2. / 5.;
/* First corner is at north pole: */
    v[12] = 0.;
    v[24] = 0.;
    v[36] = 1.;
/* The next five lie on a circle, with one on the y-axis: */
    z__ = .447213595;
/* This is 1/(2 sin^2(pi/5)) - 1 */
    rho = sqrt(1. - z__ * z__);
    for (i__ = 0; i__ <= 4; ++i__) {
	v[i__ + 13] = -rho * sin(i__ * dphi);
	v[i__ + 25] = rho * cos(i__ * dphi);
	v[i__ + 37] = z__;
    }
/* The 2nd half are simply opposite the first half: */
    for (i__ = 0; i__ <= 5; ++i__) {
	v[i__ + 18] = -v[i__ + 12];
	v[i__ + 30] = -v[i__ + 24];
	v[i__ + 42] = -v[i__ + 36];
    }
    return 0;
} /* compute_corners__ */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!     THE SUBROUTINES BELOW ARE SUBORDINATE TO THOSE ABOVE, AND  !!! */
/* !!     CAN BE SAFELY IGNORED BY THE GENERAL USER OF THE PACKAGE.  !!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!	These subroutines perform some standard linear algebra:    !!! */
/* !!	  subroutine matmul1(A,B,C)                                !!! */
/* !!	  subroutine matmul2(A,B,C)                                !!! */
/* !!	  subroutine matmul3(A,B,C)                                !!! */
/* !!	  subroutine vecmatmul1(A,b,c)	                           !!! */
/* !!	  subroutine vecmatmul2(A,b,c)	                           !!! */
/* !!                                                                !!! */
/* !!	These subroutines copy matrices in and out of storage:     !!! */
/* !!	  subroutine getmatrix(n,R,A)                              !!! */
/* !!	  subroutine putmatrix(n,R,A)                              !!! */
/* !!                                                                !!! */
/* !!     These subroutines help vector2pixel reduce the 3D sphere   !!! */
/* !!     problem to a problem on an equilateral triangle in the     !!! */
/* !!     z=1 tangent plane (an icosahedron face):                   !!! */
/* !!	  subroutine find_face(vector,R,face)                      !!! */
/* !!	  subroutine find_another_face(vector,R,face)              !!! */
/* !!	  subroutine find_corner(vector,v,corner)                  !!! */
/* !!                                                                !!! */
/* !!     These subroutines pixelize this triangle with a regular    !!! */
/* !!     triangular grid:                                           !!! */
/* !!	  subroutine find_mn(pixel,resolution,m,n)                 !!! */
/* !!	  subroutine tangentplanepixel(resolution,x,y,pix,ifail)   !!! */
/* !!	  subroutine tangentplanevector(pix,resolution,x,y)        !!! */
/* !!                                                                !!! */
/* !!     These subroutines reduce the area equalization problem to  !!! */
/* !!     one on the right triangle in the lower right corner:       !!! */
/* !!	  subroutine find_sixth(x,y,rot,flip)                      !!! */
/* !!	  subroutine rotate_and_flip(rot,flip,x,y)	           !!! */
/* !!	  subroutine adjust(x,y)                                   !!! */
/* !!	  subroutine unadjust(x,y)                                 !!! */
/* !!                                                                !!! */
/* !!     These subroutines perform the area equalization mappings   !!! */
/* !!     on this right triangle:                                    !!! */
/* !!	  subroutine adjust_sixth(x,y)                             !!! */
/* !!	  subroutine unadjust_sixth(x,y)                           !!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ****************************************************************** */
/* Subroutine */ int matmul1_(doublereal *a, doublereal *b, doublereal *c__)
{
    static doublereal d__[9]	/* was [3][3] */;
    static integer i__, j, k;
    extern /* Subroutine */ int copymatrix_(doublereal *, doublereal *);
    static doublereal sum;

/* Matrix multiplication C = AB. */
/* A, B and C are allowed to be physiclly the same. */
    /* Parameter adjustments */
    c__ -= 4;
    b -= 4;
    a -= 4;

    /* Function Body */
    sum = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    sum = 0.;
	    for (k = 1; k <= 3; ++k) {
		sum += a[i__ + k * 3] * b[k + j * 3];
	    }
	    d__[i__ + j * 3 - 4] = sum;
	}
    }
    copymatrix_(d__, &c__[4]);
    return 0;
} /* matmul1_ */

/* ****************************************************************** */
/* Subroutine */ int matmul2_(doublereal *a, doublereal *b, doublereal *c__)
{
    static doublereal d__[9]	/* was [3][3] */;
    static integer i__, j, k;
    extern /* Subroutine */ int copymatrix_(doublereal *, doublereal *);
    static doublereal sum;

/* Matrix multiplication C = AB^t */
/* A, B and C are allowed to be physically the same. */
    /* Parameter adjustments */
    c__ -= 4;
    b -= 4;
    a -= 4;

    /* Function Body */
    sum = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    sum = 0.;
	    for (k = 1; k <= 3; ++k) {
		sum += a[i__ + k * 3] * b[j + k * 3];
	    }
	    d__[i__ + j * 3 - 4] = sum;
	}
    }
    copymatrix_(d__, &c__[4]);
    return 0;
} /* matmul2_ */

/* ****************************************************************** */
/* Subroutine */ int matmul3_(doublereal *a, doublereal *b, doublereal *c__)
{
    static doublereal d__[9]	/* was [3][3] */;
    static integer i__, j, k;
    extern /* Subroutine */ int copymatrix_(doublereal *, doublereal *);
    static doublereal sum;

/* Matrix multiplication C = A^t B */
/* A, B and C are allowed to be physically the same. */
    /* Parameter adjustments */
    c__ -= 4;
    b -= 4;
    a -= 4;

    /* Function Body */
    sum = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    sum = 0.;
	    for (k = 1; k <= 3; ++k) {
		sum += a[k + i__ * 3] * b[k + j * 3];
	    }
	    d__[i__ + j * 3 - 4] = sum;
	}
    }
    copymatrix_(d__, &c__[4]);
    return 0;
} /* matmul3_ */

/* ****************************************************************** */
/* Subroutine */ int vecmatmul1_(doublereal *a, doublereal *b, doublereal *
	c__)
{
    static doublereal d__[3];
    static integer i__, j;
    extern /* Subroutine */ int copyvector_(doublereal *, doublereal *);
    static doublereal sum;

/* Matrix multiplication c = Ab */
/* b and c are allowed to be physically the same. */
    /* Parameter adjustments */
    --c__;
    --b;
    a -= 4;

    /* Function Body */
    sum = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	sum = 0.;
	for (j = 1; j <= 3; ++j) {
	    sum += a[i__ + j * 3] * b[j];
	}
	d__[i__ - 1] = sum;
    }
    copyvector_(d__, &c__[1]);
    return 0;
} /* vecmatmul1_ */

/* ****************************************************************** */
/* Subroutine */ int vecmatmul2_(doublereal *a, doublereal *b, doublereal *
	c__)
{
    static doublereal d__[3];
    static integer i__, j;
    extern /* Subroutine */ int copyvector_(doublereal *, doublereal *);
    static doublereal sum;

/* Matrix multiplication c = A^tb */
/* b and c are allowed to be physiclly the same. */
    /* Parameter adjustments */
    --c__;
    --b;
    a -= 4;

    /* Function Body */
    sum = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	sum = 0.f;
	for (j = 1; j <= 3; ++j) {
	    sum += a[j + i__ * 3] * b[j];
	}
	d__[i__ - 1] = sum;
    }
    copyvector_(d__, &c__[1]);
    return 0;
} /* vecmatmul2_ */

/* ****************************************************************** */
/* Subroutine */ int copymatrix_(doublereal *a, doublereal *b)
{
    static integer i__, j;

/* B = A */
    /* Parameter adjustments */
    b -= 4;
    a -= 4;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    b[i__ + j * 3] = a[i__ + j * 3];
	}
    }
    return 0;
} /* copymatrix_ */

/* ****************************************************************** */
/* Subroutine */ int copyvector_(doublereal *a, doublereal *b)
{
    static integer i__;

/* b = a */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	b[i__] = a[i__];
    }
    return 0;
} /* copyvector_ */

/* ****************************************************************** */
/* Subroutine */ int getmatrix_(integer *n, doublereal *r__, doublereal *a)
{
    static integer i__, j;

/* A = the nth matrix in R */
    /* Parameter adjustments */
    a -= 4;
    r__ -= 80;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    a[i__ + j * 3] = r__[*n + (i__ + j * 3) * 20];
	}
    }
    return 0;
} /* getmatrix_ */

/* ****************************************************************** */
/* Subroutine */ int putmatrix_(integer *n, doublereal *r__, doublereal *a)
{
    static integer i__, j;

/* the nth matrix in R = A */
    /* Parameter adjustments */
    a -= 4;
    r__ -= 80;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    r__[*n + (i__ + j * 3) * 20] = a[i__ + j * 3];
	}
    }
    return 0;
} /* putmatrix_ */

/* ****************************************************************** */
/* Subroutine */ int find_face__(doublereal *vector, doublereal *r__, integer 
	*face)
{
    static integer i__, n;
    static doublereal max__, dot;

/* Locates the face to which vector points. */
/* Computes the dot product with the vectors */
/* pointing to the center of each face and picks the */
/* largest one. */
/* This simple routine can be substantially accelerated */
/* by adding a bunch of if-statements, to avoid looping */
/* over more than a few faces. */
    /* Parameter adjustments */
    r__ -= 80;
    --vector;

    /* Function Body */
    max__ = -17.;
    for (n = 0; n <= 19; ++n) {
	dot = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    dot += r__[n + (i__ + 9) * 20] * vector[i__];
	}
	if (dot > max__) {
	    *face = n;
	    max__ = dot;
	}
    }
    return 0;
} /* find_face__ */

/* ****************************************************************** */
/* Subroutine */ int find_another_face__(doublereal *vector, doublereal *r__, 
	integer *face)
{
    static integer i__, n;
    static doublereal max__, dot;
    static integer facetoavoid;

/* Computes the dot product with the vectors */
/* pointing to the center of each face and picks the */
/* largest one other than face. */
/* This simple routine can be substantially accelerated */
/* by adding a bunch of if-statements, to avoid looping */
/* over more than a few faces. */
    /* Parameter adjustments */
    r__ -= 80;
    --vector;

    /* Function Body */
    facetoavoid = *face;
    max__ = -17.;
    for (n = 0; n <= 19; ++n) {
	if (n != facetoavoid) {
	    dot = 0.;
	    for (i__ = 1; i__ <= 3; ++i__) {
		dot += r__[n + (i__ + 9) * 20] * vector[i__];
	    }
	    if (dot > max__) {
		*face = n;
		max__ = dot;
	    }
	}
    }
    return 0;
} /* find_another_face__ */

/* ****************************************************************** */
/* Subroutine */ int find_corner__(doublereal *vector, doublereal *v, integer 
	*corner)
{
    static integer i__, n;
    static doublereal max__, dot;

/* Locates the corner to which vector points. */
/* Computes the dot product with the vectors */
/* pointing to each corner and picks the */
/* largest one. */
/* This simple routine can be substantially accelerated */
/* by adding a bunch of if-statements, but that's pretty */
/* pointless since it gets called so rarely. */
    /* Parameter adjustments */
    v -= 12;
    --vector;

    /* Function Body */
    max__ = -17.;
    for (n = 0; n <= 11; ++n) {
	dot = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    dot += v[n + i__ * 12] * vector[i__];
	}
	if (dot > max__) {
	    *corner = n;
	    max__ = dot;
	}
    }
    return 0;
} /* find_corner__ */

/* ****************************************************************** */
/* Subroutine */ int find_mn__(integer *pixel, integer *resolution, integer *
	m, integer *n)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer pixperedge, pix, interiorpix;

/* Computes the integer coordinates (m,n) of the pixel */
/* numbered pix on the basic triangle. */
    pix = *pixel;
    interiorpix = ((*resolution << 1) - 3) * (*resolution - 1);
    pixperedge = *resolution - 1;
    if (pix < interiorpix) {
/* The pixel lies in the interior of the triangle. */
	*m = (integer) ((sqrt(pix * 8. + 1.) - 1.) / 2. + .5 / *resolution);
/* 0.5/resolution was added to avoid problems with */
/* rounding errors for the case when n=0. */
/* As long as you don't add more than 2/m, you're OK. */
	*n = pix - *m * (*m + 1) / 2;
	*m += 2;
	++(*n);
	goto L555;
    }
    pix -= interiorpix;
    if (pix < pixperedge) {
/* The pixel lies on the bottom edge. */
	*m = (*resolution << 1) - 1;
	*n = pix + 1;
	goto L555;
    }
    pix -= pixperedge;
    if (pix < pixperedge) {
/* The pixel lies on the right edge. */
	*m = (*resolution << 1) - (pix + 2);
	*n = *m;
	goto L555;
    }
    pix -= pixperedge;
/* The pixel lies on the left edge. */
    *m = pix + 1;
    *n = 0;
L555:
    return 0;
} /* find_mn__ */

/* ****************************************************************** */
/* Subroutine */ int tangentplanepixel_(integer *resolution, doublereal *x, 
	doublereal *y, integer *pix, integer *ifail)
{
    static doublereal a, b, d__;
    static integer i__, j, k, m, n, r2;

/* Finds the hexagon in which the point (x,y) lies */
/* and computes the corresponding pixel number pix. */
/* Returns ifail=0 if (x,y) lies on the face, */
/* otherwise returns ifail=1. */
/* sqrt(3.d0)/2.d0 */
/* The edge length of the icosahedron is */
/* sqrt(9 tan^2(pi/5) - 3) when scaled so that */
/* it circumscribes the unit sphere. */
    r2 = *resolution << 1;
    a = *x * .5;
    b = *y * .866025404;
    d__ = .66158453824999996 / r2;
    i__ = (integer) (*x / d__ + r2);
    j = (integer) ((a + b) / d__ + r2);
    k = (integer) ((a - b) / d__ + r2);
    m = (r2 + r2 - j + k - 1) / 3;
    n = (i__ + k + 1 - r2) / 3;
    *pix = (m - 2) * (m - 1) / 2 + (n - 1);
    *ifail = 0;
    if (m == r2 - 1) {
/* On bottom row */
	if (n <= 0 || n >= *resolution) {
	    *ifail = 1;
	}
	goto L666;
/* Pix already correct */
    }
    if (n == m) {
/* On right edge */
	k = r2 - 1 - m;
	if (k <= 0 || k >= *resolution) {
	    *ifail = 1;
	} else {
	    *pix = (r2 - 2) * (*resolution - 1) + k - 1;
	}
	goto L666;
    }
    if (n == 0) {
/* On left edge */
	if (m <= 0 || m >= *resolution) {
	    *ifail = 1;
	} else {
	    *pix = (r2 - 1) * (*resolution - 1) + m - 1;
	}
    }
L666:
    return 0;
} /* tangentplanepixel_ */

/* ****************************************************************** */
/* Subroutine */ int tangentplanevector_(integer *pix, integer *resolution, 
	doublereal *x, doublereal *y)
{
    static integer m, n;
    extern /* Subroutine */ int find_mn__(integer *, integer *, integer *, 
	    integer *);

/* Computes the coordinates (x,y) of the pixel */
/* numbered pix on the basic triangle. */
/* 1.d0/sqrt(3.d0) */
/* sqrt(3.d0)/2.d0 */
/* The edge length of the icosahedron is */
/* sqrt(9 tan^2(pi/5) - 3) when scaled so that */
/* it circumscribes the unit sphere. */
    find_mn__(pix, resolution, &m, &n);
    *x = (n - m * .5) * 1.3231690765 / (*resolution * 2. - 1.);
    *y = (.577350269 - .866025404 / (*resolution * 2. - 1.) * m) * 
	    1.3231690765;
    return 0;
} /* tangentplanevector_ */

/* ****************************************************************** */
/* Subroutine */ int find_sixth__(doublereal *x, doublereal *y, integer *rot, 
	integer *flip)
{
    static doublereal d__;

/* Find out in which sixth of the basic triangle */
/* the point (x,y) lies, identified by the */
/* two integers rot (=0, 1 or 2) and flip = 0 or 1). */
/* rot and flip are defined such that the sixth is */
/* mapped onto the one at the bottom right by */
/* these two steps: */
/* 1. Rotate by 120 degrees anti-clockwise, rot times. */
/* 2. Flip the sign of x if flip = 1, not if flip=0. */
/* The if-statements below go through the six cases */
/* anti-clockwise, starting at the bottom right. */
/* sqrt(3d0) */
    d__ = *y * 1.73205081;
    if (*x >= 0.) {
	if (*x <= -d__) {
	    *rot = 0;
	    *flip = 0;
	} else {
	    if (*x >= d__) {
		*rot = 2;
		*flip = 1;
	    } else {
		*rot = 2;
		*flip = 0;
	    }
	}
    } else {
	if (*x >= -d__) {
	    *rot = 1;
	    *flip = 1;
	} else {
	    if (*x <= d__) {
		*rot = 1;
		*flip = 0;
	    } else {
		*rot = 0;
		*flip = 1;
	    }
	}
    }
    return 0;
} /* find_sixth__ */

/* ****************************************************************** */
/* Subroutine */ int rotate_and_flip__(integer *rot, integer *flip, 
	doublereal *x, doublereal *y)
{
    static doublereal x1, sn;

/* sqrt(3d0)/2d0 */
    if (*rot > 0) {
	if (*rot == 1) {
	    sn = .866025404;
/* Rotate 120 degrees anti-clockwise */
	} else {
	    sn = -.866025404;
/* Rotate 120 degrees anti-clockwise */
	}
	x1 = *x;
	*x = x1 * -.5 - sn * *y;
	*y = sn * x1 + *y * -.5;
    }
    if (*flip > 0) {
	*x = -(*x);
    }
    return 0;
} /* rotate_and_flip__ */

/* ****************************************************************** */
/* Subroutine */ int adjust_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int find_sixth__(doublereal *, doublereal *, 
	    integer *, integer *);
    static integer rot, flip;
    extern /* Subroutine */ int adjust_sixth__(doublereal *, doublereal *), 
	    rotate_and_flip__(integer *, integer *, doublereal *, doublereal *
	    );

/* Maps the basic triangle onto itself in such a way */
/* that pixels will have equal area when mapped onto */
/* the sphere. */
    find_sixth__(x, y, &rot, &flip);
    rotate_and_flip__(&rot, &flip, x, y);
    adjust_sixth__(x, y);
/* Now rotate & flip the sixth back into its */
/* original position: */
    if (flip == 0 && rot > 0) {
	i__1 = 3 - rot;
	rotate_and_flip__(&i__1, &flip, x, y);
    } else {
	rotate_and_flip__(&rot, &flip, x, y);
    }
    return 0;
} /* adjust_ */

/* ****************************************************************** */
/* Subroutine */ int unadjust_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int find_sixth__(doublereal *, doublereal *, 
	    integer *, integer *);
    static integer rot, flip;
    extern /* Subroutine */ int unadjust_sixth__(doublereal *, doublereal *), 
	    rotate_and_flip__(integer *, integer *, doublereal *, doublereal *
	    );

/* Performs the inverse of what adjust does. */
    find_sixth__(x, y, &rot, &flip);
    rotate_and_flip__(&rot, &flip, x, y);
    unadjust_sixth__(x, y);
/* Now rotate & flip the sixth back into its */
/* original position: */
    if (flip == 0 && rot > 0) {
	i__1 = 3 - rot;
	rotate_and_flip__(&i__1, &flip, x, y);
    } else {
	rotate_and_flip__(&rot, &flip, x, y);
    }
    return 0;
} /* unadjust_ */

/* ****************************************************************** */
/* Subroutine */ int adjust_sixth__(doublereal *x, doublereal *y)
{
    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal u, v, v2, trig, root;

/* Maps the basic right triangle (the sixth of the face that */
/* is in the lower right corner) onto itself in such a way */
/* that pixels will have equal area when mapped onto the sphere. */
/* sqrt(3.d0) */
    u = *x + 1e-14;
    v = -(*y) + 1e-14;
    v2 = v * v;
    root = sqrt(v2 * 4. + 1.);
    trig = atan((root * 1.7320508075689 - 1.7320508075689) / (root + 3.));
    *y = sqrt(trig * 2. / 1.7320508075689);
    *x = sqrt((v2 * 4. + 1.) / (u * u + 1. + v2)) * u * *y / v;
    *x *= 1.09844;
    *y *= -1.09844;
    return 0;
} /* adjust_sixth__ */

/* ****************************************************************** */
/* Subroutine */ int unadjust_sixth__(doublereal *x, doublereal *y)
{
    /* Builtin functions */
    double tan(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal u, v, v2, y2, tmp, trig;

/* Performs the inverse of what adjust_sixth does. */
/* sqrt(3.d0) */
    u = *x / 1.09844 + 1e-14;
    v = -(*y) / 1.09844 + 1e-14;
    v2 = v * v;
    trig = tan(v2 * 1.7320508075689 / 2.);
    tmp = (trig * 3. + 1.7320508075689) / (1.7320508075689 - trig);
    y2 = (tmp * tmp - 1.) / 4.;
    *y = sqrt(y2);
    tmp = v2 * (y2 * 4. + 1.) - u * u * y2;
    *x = u * *y * sqrt((y2 + 1.) / tmp);
    *y = -(*y);
    return 0;
} /* unadjust_sixth__ */

