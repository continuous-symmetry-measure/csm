#include "newton_solver.h"
#include "equation_set.h"
#include "MasterEquations.h"
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
			pn.types[i].cutoff = (size_t)floor(5 * finalState[i] + 4);
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
	if (argc == 3 && strcmp(argv[2], "-force")) {
		force = true;
	}

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

	if (!force) {
		cout << "Running Rate Equations and updating cutoffs" << endl;

		// First - run rate equations to compute avarages
 		RateEquations req(pn);
		CutoffProcessor cp;
		cp.pn = pn;
		RKSolver<double> rateSolver(cp, req);
		RateEquations::vec initial = req.createInitialConditions();		
		rateSolver.solve(params, initial);	
	
		cout << "Running Master equations with updated cutoffs " << endl;	
		MasterEquations eq(cp.pn);
		NewtonSolver<double> solver(eq, eq);
		MasterEquations::vec initialState = eq.createInitialConditions();
		solver.solve(params, initialState, 20);
	} else {
		MasterEquations eq(pn);
		NewtonSolver<double> solver(eq, eq);
		MasterEquations::vec initialState = eq.createInitialConditions();
		solver.solve(params, initialState,20);
	}
}
