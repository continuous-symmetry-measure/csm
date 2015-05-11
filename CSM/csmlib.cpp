/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>
#include <cstdio>
#include "Molecule.h"

using namespace std;

extern int main(int argc, char *argv[]); // Defined in mainRot.cpp

static csm_options process_bridge(const python_cpp_bridge &bridge);
static FILE *convert_to_file(int fd, bool *flag=NULL);

python_cpp_bridge::python_cpp_bridge()
{
	printNorm = printLocal = writeOpenu = ignoreHy = removeHy = findPerm = useMass = limitRun = babelBond = timeOnly = detectOutliers = babelTest = keepCenter = false;
	sn_max = 8;
	fdOut = -1;  // -1 means no file
}

FILE *convert_to_file(int fd, const char *mode, bool *flag=NULL)
{
#pragma warning (disable: 4996)
	FILE *f = NULL;
	if (fd != -1)
		f = fdopen(fd, mode);  // This causes a warning on Windows, but is required on Linux since _fdopen is nowhere to be found

	if (flag != NULL)
		*flag = f != NULL;

	return f;
}

csm_options process_bridge(const python_cpp_bridge &bridge)
{
	csm_options options;

	options.printNorm = bridge.printNorm;
	options.printLocal = bridge.printLocal;
	options.writeOpenu = bridge.writeOpenu;
	options.ignoreHy = bridge.ignoreHy;
	options.removeHy = bridge.removeHy;
	options.ignoreSym = bridge.ignoreSym;
	options.findPerm = bridge.findPerm;
	options.useMass = bridge.useMass;
	options.limitRun = bridge.limitRun;
	options.babelBond = bridge.babelBond;
	options.timeOnly = bridge.timeOnly;
	options.detectOutliers = bridge.detectOutliers;
	options.babelTest = bridge.babelTest;
	options.keepCenter = bridge.keepCenter;
	options.sn_max = bridge.sn_max;

	options.format = bridge.format;
	options.useFormat = options.format != "";

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

	options.inFileName = bridge.inFilename;
	options.outFileName = bridge.outFilename;
	options.logFileName = bridge.logFilename;

	options.outFile = convert_to_file(bridge.fdOut, "w");
	
	options.dir = bridge.dir;
	options.useDir = bridge.dir.size() == 3;

	options.perm = bridge.perm;
	options.useperm = bridge.perm.size() > 0;

	options.molecule = Molecule::createFromPython(bridge.molecule);

	return options;
}

extern csm_options options;
extern int mainWithOptions();

int RunCSM(python_cpp_bridge bridge)
{
	options = process_bridge(bridge);

	cout << "C++ domain entered" << endl;
	int rc = mainWithOptions();
	cout << "C++ is done" << endl;

	return 1;
}

int SayHello()
{
	cout << "Hello from C++" << endl;
	return 17;
}
