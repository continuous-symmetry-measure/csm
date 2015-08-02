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

