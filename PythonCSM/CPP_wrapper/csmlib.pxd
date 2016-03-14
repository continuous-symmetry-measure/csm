from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "csmlib.h":
    cdef cppclass python_atom:
        python_atom();
        string symbol;
        vector[int] adjacent;
        vector[double] pos;
        double mass;

cdef extern from "csmlib.h":
    cdef cppclass python_molecule:
        vector[python_atom] atoms;
        vector[vector[int]] equivalenceClasses;

cdef extern from "csmlib.h":
    cdef cppclass python_cpp_bridge:
        python_cpp_bridge()
        string opType;
        string opName;
        bool writeOpenu;

        int opOrder;
        int sn_max;
        bool detectOutliers;
        string logFilename;

        vector[double] dir;
        vector[int] perm;
        python_molecule molecule;

cdef extern from "csmlib.h":
    cdef cppclass csm_calculation_data:
        python_molecule molecule;
        vector[vector[double]] outAtoms;  # x,y,z of each atom
        vector[double] dir;
        double csm;
        double dMin;
        vector[int] perm;
        vector[double] localCSM;
        string operationType;
        int chMinOrder;
        string chMinType;
        int opOrder;

cdef extern from "csmlib.h":
    void SetCSMOptions(python_cpp_bridge options) except +;
    vector[vector[int]] GetPermuterPermutations(int size, int groupSize, bool addGroupsOfTwo);
    vector[vector[int]] GetMoleculePermutations();

    csm_calculation_data RunSinglePerm(csm_calculation_data input) except +;
    csm_calculation_data FindBestPermUsingDir (csm_calculation_data input) except +;
    csm_calculation_data FindBestPerm (csm_calculation_data input) except +;
    csm_calculation_data ComputeLocalCSM (csm_calculation_data input) except +;
    csm_calculation_data CalcRefPlane (csm_calculation_data input) except +;
    csm_calculation_data CreateSymmetricStructure (csm_calculation_data input) except +;
    int rpoly(double *coeffs, int degree, double *zeror, double *zeroi);
    void GetEigens(const double matrix[3][3], double eigenVectors[3][3], double eigenValues[3]);