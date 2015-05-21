/*
 * The interface of the CSM library.
 * The CSM library is going to be called from Python using Cython. The library's outside interface is going
 * to be located in this file, and call internal functions as necessary.
 *
 * This file is going to change continiously, until we end up with the core C++ calculations and a big Python codebase.
 *
 * By Itay Zandbank
 */

#ifndef CSMLIB_H
#define CSMLIB_H

#include <vector>
#include <string>
#include "options.h"


// A representation of one atom in Python
struct python_atom
{
	std::string symbol;
	std::vector<int> adjacent;
	std::vector<double> pos;
	double mass;

	python_atom() : mass(0.0)
	{ }
};

// Input from Python
struct python_cpp_bridge
{
	std::string opType;
	std::string opName;
	int opOrder;

	bool printNorm;
	bool printLocal;
	bool writeOpenu;
	std::string format;

	bool ignoreHy;
	bool removeHy;
	bool ignoreSym;
	bool findPerm;
	bool useMass;
	bool limitRun;
	bool babelBond;
	bool timeOnly;
	int sn_max;
	bool detectOutliers;
	bool babelTest;
	bool keepCenter;
	std::string logFilename;
	std::string inFilename;
	std::string outFilename;

	// File descriptors - -1 means no file
	int fdOut;

	// Direction Axis
	std::vector<double> dir;

	//Permutation
	std::vector<int> perm;

	// Molecule
	std::vector<python_atom> molecule;

	python_cpp_bridge();
};

struct csm_output
{
	// Molecule
	std::vector<python_atom> molecule;
	double norm;
	int numGroups;

	// Results from mainRot
	std::vector<std::vector<double>> outAtoms; // x,y,z of each atom
	double csm; // The actual CSM score
	std::vector<double> dir; 
	double dMin;
	std::vector<double> localCSM;
	int chMinOrder;
	std::vector<int> perm;
};


#ifdef __cplusplus
extern "C"
{
#endif
	// Runs the entire CSM application
	// int RunCSM(const std::vector<std::string> args);
	csm_output RunCSM(python_cpp_bridge options);
#ifdef __cplusplus
}
#endif
#endif