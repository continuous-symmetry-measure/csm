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
#include "Molecule.h"


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

struct python_molecule
{
	std::vector<python_atom> atoms;   // The atoms
	std::vector<std::vector<int> > equivalenceClasses;  // Equivalence classes
};

// Input from Python
struct python_cpp_bridge
{
	std::string opType;
	std::string opName;
	int opOrder;

	bool writeOpenu;

	int sn_max;
	bool detectOutliers;
	std::string logFilename;

	// Direction Axis
	std::vector<double> dir;

	//Permutation
	std::vector<int> perm;

	python_molecule molecule;

	python_cpp_bridge();
};

struct csm_output
{
	// Molecule
	python_molecule molecule;

	// Results from mainRot
	std::vector<std::vector<double> > outAtoms; // x,y,z of each atom
	double csm; // The actual CSM score
	std::vector<double> dir; 
	double dMin;
	std::vector<double> localCSM;
	int chMinOrder;
	std::string chMinType;
	std::vector<int> perm;
};

struct csm_calculation_data
{
	python_molecule molecule;
	std::vector<std::vector<double> > outAtoms; // x,y,z of each atom
	std::vector<double> dir;
	double csm;
	double dMin;
	std::vector<int> perm;
	std::vector<double> localCSM;
	std::string operationType;
	int chMinOrder;
	std::string chMinType;
	int opOrder;
};

struct cpp_calculation_data
{
	Molecule *molecule;
	double **outAtoms;
	double *dir;
	double csm;
	double dMin;
	int *perm;
	double *localCSM;
	OperationType operationType;
	int chMinOrder;
	OperationType chMinType;
	int opOrder;
	

	cpp_calculation_data(const csm_calculation_data &python);
	~cpp_calculation_data();

	csm_calculation_data get_csm_data();
};

// Sets the CSM options for all future function calls
void SetCSMOptions(python_cpp_bridge options);

csm_calculation_data RunSinglePerm(csm_calculation_data input);
csm_calculation_data FindBestPermUsingDir (csm_calculation_data input);
csm_calculation_data FindBestPerm (csm_calculation_data input);
csm_calculation_data ComputeLocalCSM (csm_calculation_data input);

std::vector< std::vector<int> > GetPermuterPermutations(int size, int groupSize, bool addGroupsOfTwo);
std::vector< std::vector<int> > GetMoleculePermutations();

csm_calculation_data CalcRefPlane (csm_calculation_data input);
csm_calculation_data CreateSymmetricStructure (csm_calculation_data input);

extern "C"
{
	int rpoly(double *op, int degree, double *zeror, double *zeroi);
}
#endif
