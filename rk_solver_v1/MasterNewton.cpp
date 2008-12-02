#include "rk_solver.h"
#include "equation_set.h"
#include "MasterEquations.h"
#include "math.h"
#include <iostream>

using namespace std;
int main(int argc, char *argv[]) { 
	ifstream input;
	input.open(argv[1]);
	ChemicalNetwork *cn = ChemicalNetwork::parseChemicalNetwork(input);
	input.close();
	if (cn == NULL) return 1;
	
	cn->performNewtonMethod(5e-7);
}
