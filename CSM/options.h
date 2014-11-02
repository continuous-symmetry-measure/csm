/*
 * The CSM runtime options
 *
 * Extracted from mainRot.cpp during the 2014 reorganization
 * Created by Itay Zandbank
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

namespace csm_options
{
	typedef enum {
		CN,
		SN,
		CS,
		CI,
		CH
	} OperationType;

	extern char opName[100];
	extern bool printNorm;
	extern bool printLocal;
	extern bool writeOpenu;
	extern std::string format;

	extern bool ignoreHy;
	extern bool removeHy;
	extern bool ignoreSym;
	extern bool useFormat;
	extern OperationType type;
	extern int opOrder;
	extern bool useperm;
	extern bool useDir;
	extern bool findPerm;
	extern bool useMass;
	extern bool limitRun;
	extern bool babelBond;
	extern bool timeOnly;
	extern int sn_max;
	extern bool detectOutliers;
	extern double A;
	extern bool babelTest;
	extern bool keepCenter;
	extern std::string logFile;

	// file pointers
	extern FILE* inFile;
	extern FILE* outFile;
	extern FILE* permfile;
	extern FILE* dirfile;
	extern char *inFileName;
	extern char *outFileName;

	void usage(const std::string op);
	void parseInput(int argc, char *argv[]);
}
#endif