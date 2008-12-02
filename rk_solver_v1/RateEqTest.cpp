#include "rk_solver.h"
#include "equation_set.h"
#include "math.h"
#include <iostream>

using namespace std;

double d;
#define T	(10)
#define NU	(1e12)
#define s	(5e13)	// site density per micrometer
#define k	(8.617343e-5)

#define	S	(3.14 * d * d * s)
#define F	(1e-11 * S)
#define W	(NU * exp(-32e-3 / k / T))
#define A	((NU * exp(-22e-3 / k / T)) / S)

class RateEqSet : public EquationSet<double>{
	virtual vec compute_derivative(const double time, const vec& state) {
		vec res = state.copy();
		res[0] = F - W * state[0] - 2 * A * state[0] * state[0];
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
		file << time << " " << dt << " H: " << state[0] << endl;
	}

	/** 
	* Announce that the solving is complete
	*/
	virtual void solutionComplete(const res_vec& finalState) { 
		file << "Finished Solving" << endl;
	}

	virtual ~ResultProcessor() {}
};

int main(int argc, char *argv[]) { 
	
	cin >> d;	
	d /= 10000;

	cout << "S: " << S << endl;
	cout << "F: " << F << endl;
	cout << "W: " << W << endl;	
	cout << "A: " << A << endl;

	ResultProcessor processor_rate("rate.out");
	RateEqSet rateEqs;
	RKSolver<double> solver_rate(processor_rate, rateEqs);
	rk_params params;
	RateEqSet::vec initialStateRate(1);
	initialStateRate[0] = 0.0;
	
	params.final_time = 1000000;
	params.initial_time = 0;
	params.max_error = 1e-7;
	params.limit = 5e-6;
	params.initialDelta = 1e-10;

	ifstream is;
	is.open("rk_params");	
	if (is.good()) {
		is >> params.initial_time >> params.final_time >> params.max_error >> params.limit >> params.initialDelta;
	}
	is.close();


	solver_rate.solve(params, initialStateRate);


}
