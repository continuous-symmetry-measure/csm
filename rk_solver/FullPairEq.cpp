#include "rk_solver.h"
#include "equation_set.h"
#include "FullPairEquations.h"
#include "RateEquations.h"
#include "math.h"
#include <iostream>

using namespace std;

class CutoffProcessor : public RKResultProcessor<double> {
public:
	parsed_network pn;

	virtual void solutionComplete(const res_vec& finalState) { 
		cout << endl;
		for (size_t i = 0; i < finalState.size(); i++) {
			if (floor(5 * finalState[i] + 5 + log(finalState.max()) / log(10)) > 2) {
				pn.types[i].cutoff = (size_t)floor(5 * finalState[i] + 5 + log(finalState.max()) / log(10));
			} else {
				pn.types[i].cutoff = 2;
			}
			cout << "Updated Cutoff for " << pn.types[i].name << " is " << pn.types[i].cutoff << endl;
		}
	}
};

int main(int argc, char *argv[]) {
	ifstream input;
	input.open(argv[1]);
	parsed_network pn = ChemicalNetwork::parseChemicalNetwork(input);
	input.close();

	bool force = false;
	if (argc == 3 && strcmp(argv[2], "-force") == 0) {
		force = true;
	}

	rk_params params;
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

	if (!force) {
		cout << "Running Rate Equations and updating cutoffs" << endl;

		// First - run rate equations to compute avarages
 		RateEquations req(pn);
		CutoffProcessor cp;
		cp.pn = pn;
		RKSolver<double> rateSolver(cp, req);
		RateEquations::vec initial = req.createInitialConditions();		
		rateSolver.solve(params, initial);	
	
		cout << "Running Pair equations with updated cutoffs " << endl;	
		FullPairEquations eq(cp.pn);
		RKSolver<Matrix<double> > solver(eq, eq);
		FullPairEquations::vec initialState = eq.createInitialConditions();
		solver.solve(params, initialState);
	} else {
		FullPairEquations eq(pn);
		FullPairEquations::vec initialState = eq.createInitialConditions();
		solver.solve(params, initialState);
	}
}