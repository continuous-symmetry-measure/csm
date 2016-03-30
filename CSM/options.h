/*
 * The CSM runtime options
 *
 * Extracted from mainRot.cpp during the 2014 reorganization
 * Created by Itay Zandbank
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <vector>

typedef enum {
	CN,
	SN,
	CS,
	CI,
	CH
} OperationType;

class Molecule;  // Forward reference

struct csm_options
{
	std::string opName;
	bool writeOpenu;

	OperationType type;
	int opOrder;
	int sn_max;
	bool detectOutliers;
	std::string logFileName;

	std::vector<double> dir;
	std::vector<int> perm;

	Molecule *molecule;  // The molecule to be used

	csm_options();
};

extern csm_options options;

#endif