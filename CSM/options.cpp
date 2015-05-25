/*
* The CSM runtime options
*
* Extracted from mainRot.cpp during the 2014 reorganization
* Created by Itay Zandbank
*/

#include "options.h"

using namespace std;

csm_options::csm_options() 
{ 
	printNorm = printLocal = writeOpenu = ignoreHy = removeHy = useFormat = useperm = useDir = findPerm = useMass = limitRun = babelBond = timeOnly = detectOutliers = babelTest = keepCenter = false;
	outFile = NULL;
	dir.clear();
	perm.clear();
	molecule = NULL;
}