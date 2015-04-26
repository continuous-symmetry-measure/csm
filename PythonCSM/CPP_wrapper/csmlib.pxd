from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "csmlib.h":
    cdef cppclass python_cpp_bridge:
        python_cpp_bridge()
        string opType;
        string opName;
        bool printNorm;
        bool printLocal;
        bool writeOpenu;
        string format;

        bool ignoreHy;
        bool removeHy;
        bool ignoreSym;
        bool useFormat;

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
        string logFilename;

        string inFilename;
        string outFilename;

        int fdIn, fdOut, fdPerm, fdDir;

cdef extern from "csmlib.h":
    int SayHello();
    int RunCSM(python_cpp_bridge options);

