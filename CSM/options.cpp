/*
* The CSM runtime options
*
* Extracted from mainRot.cpp during the 2014 reorganization
* Created by Itay Zandbank
*/

#include "options.h"
#include "logging.h"
#include <cstring>
#include <stdio.h>

using namespace std;

csm_options::csm_options() 
{ 
	printLocal = writeOpenu = useperm = useDir = findPerm = limitRun = timeOnly = detectOutliers = babelTest = displayPerms= false;
	dir.clear();
	perm.clear();
	molecule = NULL;
}