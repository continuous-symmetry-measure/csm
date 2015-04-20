from libcpp.string cimport string

cdef extern from "options.h":
    cdef cppclass csm_options:
        csm_options()
       	string opName;
        bool printNorm;
        bool printLocal;
        bool writeOpenu;
        string format;

        bool ignoreHy;
        bool removeHy;
        bool ignoreSym;
        bool useFormat;
        #OperationType type;
        int opOrder;
        bool useperm;
        bool useDir;
        bool findPerm;
        bool useMass;
        bool limitRun;
        bool babelBond;
        bool timeOnly;
        int sn_max;
        bool detectOutliers;
        double A;
        bool babelTest;
        bool keepCenter;
        string logFile;

        int inFile;
        int outFile;
        int permfile;
        int dirfile;

        string inFileName;
        string outFileName;

cdef extern from "csmlib.h":
    int SayHello();
    int RunCSM(csm_options options);

