#include "rk_solver.h"
#include "equation_set.h"
#include "PairEquations.h"
#include "math.h"
#include <iostream>

using namespace std;
int main(int argc, char *argv[]) {
	ifstream input;
	input.open(argv[1]);
	PairEquations eq(ChemicalNetwork::parseChemicalNetwork(input));
	input.close();
	RKSolver<Matrix<double> > solver(eq, eq);
	rk_params params;

	params.final_time = 3e10;
	params.initial_time = 0;
	params.max_error = 1e-8;
	params.limit = 1e-8;
	params.initialDelta = 1e-10;

	ifstream is;
	is.open("rk_params");
	if (is.good()) {
		is >> params.initial_time >> params.final_time >> params.max_error >> params.limit >> params.initialDelta;
	}
	is.close();

	PairEquations::vec initialState = eq.createInitialConditions();
	solver.solve(params, initialState);
}
