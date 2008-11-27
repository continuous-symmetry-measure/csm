#include "rk_solver.h"
#include "equation_set.h"
#include "math.h"
#include <iostream>

using namespace std;

double S;
#define F 1e-8

#define	F1 (F * S)
#define F2 (F1 * 0.01)
#define F3 (F1 * 0.01)
#define A1 (10.0 / S)
#define A2 (1.0 / S)
#define A3 (1e-1 / S)
#define W1 (1e-3)
#define W2 (1e-3)
#define W3 (1e-5)

class RateEqSet : public EquationSet<double>{
	virtual vec compute_derivative(const double time, const vec& state) {
		vec res = state.copy();
		res[0] = F1 - W1 * state[0] - (A1 + A2) * state[0] * state[1] - (A1 + A3) * state[0] * state[2];
		res[1] = F2 - W2 * state[1] - (A1 + A2) * state[0] * state[1];
		res[2] = F3 - W3 * state[2] - (A1 + A3) * state[0] * state[2];
		return res;
	}	
};


class ResultProcessor : public RKResultProcessor<double> {
private:
	ofstream file;
	int stepNum;

public:
	ResultProcessor(char *filename) {
		file.open(filename);
	}

	/**
	* Announce that the solving has started
	*/
	virtual void solvingStarted(const rk_params &params, const res_vec& initalState) { 
		file << "Solving Started" << endl;
	}

	/** 
	* A single complete RK step has been performed
	* @param time	The time
	* @param dt	The chosen delta
	* @param state The state after the step
	*/
	virtual void stepPerformed(double time, double dt, const res_vec& state, const res_vec& prevState) { 
		file << time << " " << dt << " X1: " << state[0] << " X2: " << state[1] << " X3: " << state[2] << endl;
	}

	/** 
	* Announce that the solving is complete
	*/
	virtual void solutionComplete(const res_vec& finalState) { 
		file << "Finished Solving" << endl;
		cout << " " << " X1: " << finalState[0] << " X2: " << finalState[1] << " X3: " << finalState[2] << endl;
		
	}

	virtual ~ResultProcessor() {}
};

int main(int argc, char *argv[]) { 
	
	S = atof(argv[1]);

	ResultProcessor processor_rate("rate.out");
	RateEqSet rateEqs;
	RKSolver<double> solver_rate(processor_rate, rateEqs);
	rk_params params;

	// Prepare network params file
	ofstream os("rate.net");
	
	os << "Begin" << endl
		<< "Species\tFlux\tDiffusion\tSweepRate\tCutoff" << endl 
		<< "X1\t" << F1 << "\t" << W1 << "\t" << A1 << "\t20" << endl
		<< "X2\t" << F2 << "\t" << W2 << "\t" << A2 << "\t20" << endl
		<< "X3\t" << F3 << "\t" << W3 << "\t" << A3 << "\t20" << endl
		<< "X4\t" << 0 << "\t" << 0 << "\t" << 0 << "\t20" << endl
		<< "X5\t" << 0 << "\t" << 0 << "\t" << 0 << "\t20" << endl
		<< endl
		<< "Interaction" << endl
		<< "X1 X2 => X4" << endl	
		<< "X1 X3 => X5" << endl	
		<< "End" << endl;
	os.close();
	RateEqSet::vec initialStateRate(3);
	initialStateRate = 0.0;
	
	params.final_time = 1e9;
	params.initial_time = 0;
	params.max_error = 1e-11;
	params.limit = 5e-7;
	params.initialDelta = 1e-10;

	ifstream is;
	is.open("rk_params");	
	if (is.good()) {
		is >> params.initial_time >> params.final_time >> params.max_error >> params.limit >> params.initialDelta;
	}
	is.close();
	solver_rate.solve(params, initialStateRate);


}
