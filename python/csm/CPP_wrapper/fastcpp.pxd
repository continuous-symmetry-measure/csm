
cdef extern from "FastCPPUtils.h":
    int rpoly(double *coeffs, int degree, double *zeror, double *zeroi);
    void GetEigens(const double matrix[3][3], double eigenVectors[3][3], double eigenValues[3]);