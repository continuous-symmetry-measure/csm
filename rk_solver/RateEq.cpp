#include "rk_solver.h"
#include "ChemicalNetwork.h"
#include "RateEquations.h"
#include "math.h"
#include <iostream>

using namespace std;
int main(int argc, char *argv[]) { 
	ifstream input;
	input.open(argv[1]);
	RateEquations eq(ChemicalNetwork::parseChemicalNetwork(input));
	input.close();
	if (cn == NULL) return 1;	RKSolver<double> solver(eq, eq);
	rk_params params;

	ChemicalNetwork::vec initialState = eq.createInitialConditions();	params.final_time = 3e10;
	params.initial_time = 0;
	params.max_error = 1e-8;
	params.limit = 1e-6;
	params.initialDelta = 1e-10;

	ifstream is;
	is.open("rk_params");	
	if (is.good()) {
		is >> params.initial_time >> params.final_time >> params.max_error >> params.limit >> params.initialDelta;
	}
	is.close();

	solver.solve(params, initialState);
}

