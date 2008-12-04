#include "rk_solver.h"
#include "equation_set.h"
#include "MasterEquations.h"
#include "math.h"
#include <iostream>

using namespace std;
int main(int argc, char *argv[]) {
	ifstream input;
	input.open(argv[1]);
	MasterEquations eq(ChemicalNetwork::parseChemicalNetwork(input));
	input.close();

	eq.performNewtonMethod(5e-7);
}
