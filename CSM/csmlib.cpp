/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>
#include "Molecule.h"
#include "calculations.h"

using namespace std;

extern int main(int argc, char *argv[]); // Defined in mainRot.cpp

static csm_options process_bridge(const python_cpp_bridge &bridge);

python_cpp_bridge::python_cpp_bridge()
{
	printLocal = writeOpenu = findPerm = limitRun = timeOnly = detectOutliers = babelTest = displayPerms = false;
	sn_max = 8;
}

csm_options process_bridge(const python_cpp_bridge &bridge)
{
	csm_options options;

	options.printLocal = bridge.printLocal;
	options.writeOpenu = bridge.writeOpenu;
	options.findPerm = bridge.findPerm;
	options.limitRun = bridge.limitRun;
	options.timeOnly = bridge.timeOnly;
	options.detectOutliers = bridge.detectOutliers;
	options.babelTest = bridge.babelTest;
	options.sn_max = bridge.sn_max;
	options.displayPerms = bridge.displayPerms;

	if (bridge.opType == "CS")
		options.type = CS;
	else if (bridge.opType == "CH")
		options.type = CH;
	else if (bridge.opType == "CN")
		options.type = CN;
	else if (bridge.opType == "SN")
		options.type = SN;
	else if (bridge.opType == "CI")
		options.type = CI;
	options.opName = bridge.opName;
	options.opOrder = bridge.opOrder;

	options.logFileName = bridge.logFilename;
	
	options.dir = bridge.dir;
	options.useDir = bridge.dir.size() == 3;

	options.perm = bridge.perm;
	options.useperm = bridge.perm.size() > 0;

	options.molecule = Molecule::createFromPython(bridge.molecule);

	return options;
}

cpp_calculation_data::cpp_calculation_data(const csm_calculation_data &python)
{
	int i;
	molecule = Molecule::createFromPython(python.molecule);
	int size = molecule->size();
	
	outAtoms = (double **)malloc(size * sizeof(double*));
	for (i=0; i<size; i++)
	{
		outAtoms[i] = (double *)malloc(DIM * sizeof(double));
	}
	
	dir = (double *) malloc(DIM * sizeof(double));
	for (i = 0; i<python.dir.size(); i++)
	{
		dir[i] = python.dir[i];
	}
	
	csm = 0;
	dMin = 0;
	
	perm = (int *)malloc(size * sizeof(int));
	int perm_size = python.perm.size();
	for (i=0; i<perm_size; i++)
	{
		perm[i] = python.perm[i];
	}
	
	localCSM = (double *)malloc(sizeof(double) * size);
		
	if (python.operationType == "CS")
		operationType = CS;
	else if (python.operationType == "CH")
		operationType = CH;
	else if (python.operationType == "CN")
		operationType = CN;
	else if (python.operationType == "SN")
		operationType = SN;
	else if (python.operationType == "CI")
		operationType = CI;
}

cpp_calculation_data::~cpp_calculation_data()
{
	int i, size = molecule->size();
	delete molecule;

	for (i=0; i<size; i++)
	{
		free(outAtoms[i]);
	}
	free(outAtoms);
	free(dir);
	free(perm);
	free(localCSM);
}

csm_calculation_data cpp_calculation_data::get_csm_data()
{
	csm_calculation_data python;

	python.molecule.atoms.clear();
	python.outAtoms.clear();
	python.localCSM.clear();
	python.perm.clear();

	int size = molecule->size();

	for (int i = 0; i < size; i++)
	{
		// Molecule
		python_atom atom;
		atom.symbol = molecule->symbol(i);
		atom.mass = molecule->mass(i);
		for (int j = 0; j < 3; j++)
			atom.pos.push_back(molecule->pos()[i][j]);
		for (int j = 0; j < molecule->valency(i); j++)
			atom.adjacent.push_back(molecule->adjacent(i, j));
		python.molecule.atoms.push_back(atom);

		// outAtoms
		std::vector<double> outAtom;
		for (int j = 0; j < 3; j++)
			outAtom.push_back(outAtoms[i][j]);
		python.outAtoms.push_back(outAtom);

		// localCSM
		if (localCSM)
			python.localCSM.push_back(localCSM[i]);

		// perm
		python.perm.push_back(perm[i]);
	}

	// Molecule equivalencyClasses
	python.molecule.equivalenceClasses.clear();
	int *group = new int[size]; // No group is larger than the molecule - this is enough
	for (int i = 1; i <= molecule->groupNum(); i++)
	{
		int groupSize = molecule->getGroup(i, group);
		python.molecule.equivalenceClasses.push_back(vector<int>(group, group + groupSize));
	}
	delete[] group;


	// Other values
	python.csm = csm;

	python.dir.clear();
	for (int i = 0; i < 3; i++)
		python.dir.push_back(dir[i]);

	python.dMin = dMin;

	return python;
}
	

extern csm_options options;
extern csm_output results;
extern int mainWithOptions();

csm_output RunCSM()
{
	mainWithOptions();  // Fills output

	return results;
}

void SetCSMOptions(python_cpp_bridge bridge)
{
	options = process_bridge(bridge);
}

double TotalNumberOfPermutations()
{
	Molecule *m = options.molecule;

	return totalNumPermutations(m);  // Notice the small t - this is the original function
}

csm_calculation_data RunSinglePerm(csm_calculation_data input)
{
	
	cpp_calculation_data cpp_input(input);
	runSinglePerm(cpp_input.molecule, cpp_input.outAtoms, cpp_input.perm, &cpp_input.csm, cpp_input.dir, &cpp_input.dMin, cpp_input.operationType);

	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}
