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
        bool printLocal;
        bool writeOpenu;

        int opOrder;
        bool findPerm;
        bool limitRun;
        bool timeOnly;
        int sn_max;
        bool detectOutliers;
        bool babelTest;
        bool displayPerms;
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

cdef extern from "csmlib.h":
    cdef cppclass csm_output:
        python_molecule molecule;

        vector[vector[double]] outAtoms;
        double csm;
        vector[double] dir;
        double dMin;
        vector[double] localCSM;
        int chMinOrder;
        string chMinType;
        vector[int] perm;

cdef extern from "csmlib.h":
    void SetCSMOptions(python_cpp_bridge options) except +;
    double TotalNumberOfPermutations();
    csm_output RunCSM() except +;
    csm_calculation_data RunSinglePerm(csm_calculation_data input) except +;

