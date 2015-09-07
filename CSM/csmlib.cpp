/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>
#include "Molecule.h"
#include "calculations.h"
#include "logging.h"
#include "permuter.h"
#include <sstream>
#include <iomanip>
#include "groupPermuter.h"

using namespace std;

static csm_options process_bridge(const python_cpp_bridge &bridge);

python_cpp_bridge::python_cpp_bridge()
{
	writeOpenu = detectOutliers = false;
	sn_max = 8;
}

csm_options process_bridge(const python_cpp_bridge &bridge)
{
	csm_options options;

	options.writeOpenu = bridge.writeOpenu;
	options.detectOutliers = bridge.detectOutliers;
	options.sn_max = bridge.sn_max;

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

	options.perm = bridge.perm;

	options.molecule = Molecule::createFromPython(bridge.molecule);

	return options;
}

cpp_calculation_data::cpp_calculation_data(const csm_calculation_data &python)
{
	int i,j;
	molecule = Molecule::createFromPython(python.molecule);
	int size = molecule->size();
	
	outAtoms = (double **)malloc(size * sizeof(double*));
	for (i=0; i<size; i++)
	{
		outAtoms[i] = (double *)malloc(DIM * sizeof(double));
		if (python.outAtoms.size()>0)
			for (j=0; j<DIM; j++)
				outAtoms[i][j] = python.outAtoms[i][j];
		else
			for (j=0; j<DIM; j++)
				outAtoms[i][j] = 0.0;
	}
	
	dir = (double *) malloc(DIM * sizeof(double));
	for (i = 0; i<python.dir.size(); i++)
	{
		dir[i] = python.dir[i];
	}
	
	csm = python.csm;
	dMin = python.dMin;
	
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

	chMinOrder = python.chMinOrder;

	if (python.chMinType == "CS")
		chMinType = CS;
	else if (python.chMinType == "CH")
		chMinType = CH;
	else if (python.chMinType == "CN")
		chMinType = CN;
	else if (python.chMinType == "SN")
		chMinType = SN;
	else if (python.chMinType == "CI")
		chMinType = CI;
	
	options.opOrder = python.opOrder; 
	
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
		for (int j = 0; j < DIM; j++)
		{
			outAtom.push_back(outAtoms[i][j]);
		}	
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
	
	python.chMinOrder = chMinOrder;
	
	switch (chMinType) {
		case CN: python.chMinType = "CN"; break;
		case SN: python.chMinType = "SN"; break;
		case CS: python.chMinType = "CS"; break;
		case CI: python.chMinType = "CI"; break;
		case CH: python.chMinType = "CH"; break;
	}

	switch (operationType) {
		case CN: python.operationType = "CN"; break;
		case SN: python.operationType = "SN"; break;
		case CS: python.operationType = "CS"; break;
		case CI: python.operationType = "CI"; break;
		case CH: python.operationType = "CH"; break;
	}
	python.opOrder = options.opOrder;

	return python;
}
	

extern csm_options options;
extern csm_output results;

void SetCSMOptions(python_cpp_bridge bridge)
{
	options = process_bridge(bridge);
	
	init_logging();
	LOG(info) << "C++ CSM starting up";

	if (options.logFileName != "")
		set_file_logging(options.logFileName);
}

csm_calculation_data RunSinglePerm(csm_calculation_data input)
{
	
	cpp_calculation_data cpp_input(input);
	runSinglePerm(cpp_input.molecule, cpp_input.outAtoms, cpp_input.perm, &cpp_input.csm, cpp_input.dir, &cpp_input.dMin, cpp_input.operationType);

	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}

csm_calculation_data FindBestPermUsingDir (csm_calculation_data input)
{
	cpp_calculation_data cpp_input(input);
	findBestPermUsingDir(cpp_input.molecule, cpp_input.outAtoms, cpp_input.perm, &cpp_input.csm, cpp_input.dir, &cpp_input.dMin, cpp_input.operationType);
	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}

csm_calculation_data FindBestPerm (csm_calculation_data input)
{
	cpp_calculation_data cpp_input(input);
	findBestPerm(cpp_input.molecule, cpp_input.outAtoms, cpp_input.perm, &cpp_input.csm, cpp_input.dir, &cpp_input.dMin, cpp_input.operationType);
	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}

csm_calculation_data ComputeLocalCSM (csm_calculation_data input)
{
	cpp_calculation_data cpp_input(input);
	computeLocalCSM(cpp_input.molecule, cpp_input.localCSM, cpp_input.perm, cpp_input.dir, 
		cpp_input.operationType != CH ? cpp_input.operationType : cpp_input.chMinType);
	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}


std::vector< std::vector<int> > GetPermuterPermutations(int size, int groupSize, bool addGroupsOfTwo)
{
	std::vector< std::vector<int> > perms;

	Permuter p(size, groupSize, addGroupsOfTwo);
	while (p.next())
	{
		vector<int> perm;
		for (int i = 0; i < size; i++)
			perm.push_back(p[i]);
		perms.push_back(perm);
	}

	return perms;
}

std::vector< std::vector<int> > GetMoleculePermutations()
{
	std::vector< std::vector<int> > perms;

	Molecule *m = options.molecule;

	int *groupSizes;
	int addGroupsOfTwo;
	std::vector<int> idxToPos;

	groupSizes = new int[m->groupNum()];
	for (int i = 1; i <= m->groupNum(); i++){
		groupSizes[i - 1] = m->getGroupSize(i);
	}

	// build idxToPos
	idxToPos.resize(m->size());
	int pos = 0;
	// build idxToPos
	for (int j = 1; j <= m->groupNum(); j++){
		for (int i = 0; i< m->size(); i++){
			if (m->similar(i) == j){
				idxToPos[pos++] = i;
			}
		}
	}

	// create permuter
	if (options.type == SN && options.opOrder > 2) {
		addGroupsOfTwo = 1;
	}
	else {
		addGroupsOfTwo = 0;
	}

	GroupPermuter gp(m->groupNum(), groupSizes, m->size(), options.opOrder, addGroupsOfTwo);

	// calculate csm for each valid permutation & remember minimal (in optimalAntimer)
	int count = 0;
	while (gp.next())
	{
		vector<int> perm(m->size());
		for (int i = 0; i < m->size(); i++)
			perm[i] = idxToPos[gp[i]];
		perms.push_back(perm);
	}

	delete[] groupSizes;

	return perms;
}


csm_calculation_data CalcRefPlane (csm_calculation_data input) {
	cpp_calculation_data cpp_input(input);
	double result_csm = calcRefPlane(cpp_input.molecule, cpp_input.perm, cpp_input.dir, 
		cpp_input.operationType);
	csm_calculation_data output = cpp_input.get_csm_data();
	output.csm = result_csm;
	return output;
}

csm_calculation_data CreateSymmetricStructure (csm_calculation_data input) {
	cpp_calculation_data cpp_input(input);
	createSymmetricStructure(cpp_input.molecule, cpp_input.outAtoms, cpp_input.perm, cpp_input.dir, 
		cpp_input.operationType, cpp_input.dMin);
	csm_calculation_data output = cpp_input.get_csm_data();
	return output;
}