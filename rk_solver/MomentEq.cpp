#include "rk_solver.h"
#include "equation_set.h"
#include "math.h"
#include <iostream>

using namespace std;
double d;

#define T	(15)
#define NU	(1e12)
#define s	(5e13)	// site density per cm^2
#define k	(8.617343e-5)

#define	S	(3.1415926535897931 * d * d * s)
#define F1	(5e-10 * S)
#define F2	(0.02 * F1)
#define F3	(0)

#define W1	(NU * exp(-52e-3 / k / T))
#define W2	(NU * exp(-54e-3 / k / T))
#define W3	(NU * exp(-54e-3 / k / T))

#define A1	((NU * exp(-44e-3 / k / T)) / S)
#define A2	((NU * exp(-47e-3 / k / T)) / S)
#define A3	((NU * exp(-47e-3 / k / T)) / S)

class MomentEqSet : public EquationSet<double>{
public:
	virtual ~MomentEqSet() {}
	// State:
	// [ <N1> <N2> <N3> <N1^2> <N2^2> <N1N2> <N1N3> ] 
	virtual vec compute_derivative(const double time, const vec& state) {
		vec res = state.copy();
		res[0] = F1 - W1 * state[0] - 2 * A1 * (state[3] - state[0]) - 
			(A1 + A2) * state[5] - (A1 + A3) * state[6];

		res[1] = F2 - W2 * state[1] - 2 * A2 * (state[4] - state[1]) -
			(A1 + A2) * state[5];
		
		res[2] = F3 - W3 * state[2] - (A1 + A3) * state[6] + 
			(A1 + A2) * state[5];
		
		res[3] = F1 + (2 * F1 + W1 + 4 * A1) * state[0] -
			(2 * W1 + 4 * A1) * state[3] - 
			(A1 + A2) * state[5] - (A1 + A3) * state[6];

		res[4] = F2 + (2 * F2 + W2 + 4 * A2) * state[1] -
			(2 * W2 + 4 * A2) * state[4] - 
			(A1 + A2) * state[5];

		res[5] = F1 * state[1] + F2 * state[0] - 
			(W1 + W2 + A1 + A2) * state[5];

		res[6] = F1 * state[2] + F3 * state[0] - 
			(W1 + W3 + A1 + A3) * state[6];
			
		return res;
	}	};


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
		//file << "Solving Started" << endl;
	}

	/** 
	* A single complete RK step has been performed
	* @param time	The time
	* @param dt	The chosen delta
	* @param state The state after the step
	*/
	virtual void stepPerformed(double time, double dt, const res_vec& state, const res_vec& prevState) { 
		//file << "t: " << time 
		//	<< " dt: " << dt 
		//	<< " H: " << state[0]  
		//	<< " O: " << state[1] << endl ;
		file << time << " " << state[0] << " " << state[1] <<endl;
	}

	/** 
	* Announce that the solving is complete
	*/
	virtual void solutionComplete(const res_vec& finalState) { 
		//file << "Finished Solving" << endl;
		cout << (d * 10000) << " " << finalState[0] << " " << finalState[1] << " " << finalState[2] << endl;

	// Prepare network params file
	ofstream os("params.net");
	
	os << "Begin" << endl
		<< "Species\tFlux\tDiffusion\tSweepRate\tCutoff" << endl 
		<< "H\t" << F1 << "\t" << W1 << "\t" << A1 << "\t" << (int)(finalState[0] * 5 + 4) << endl
		<< "O\t" << F2 << "\t" << W2 << "\t" << A2 << "\t" << (int)(finalState[1] * 5 + 4) << endl
		<< "OH\t" << F3 << "\t" << W3 << "\t" << A3 << "\t" << (int)(finalState[2] * 5 + 4) << endl
		<< endl
		<< "Interaction" << endl
		<< "H H => H2" << endl	
		<< "H O => OH" << endl
		<< "H OH => H2O" << endl
		<< "O O => O2" << endl
		<< "End" << endl;
	os.close();
	}

	virtual ~ResultProcessor() {}
};

int main(int argc, char *argv[]) { 
	sscanf(argv[1], "%lf", &d);
	d /= 10000;

	cout << "S: " << S << endl;
	cout << "F1: " << F1 << endl;
	cout << "W1: " << W1 << endl;	
	cout << "A1: " << A1 << endl;

	cout << "F2: " << F2 << endl;
	cout << "W2: " << W2 << endl;	
	cout << "A2: " << A2 << endl;

	cout << "F3: " << F3 << endl;
	cout << "W3: " << W3 << endl;	
	cout << "A3: " << A3 << endl;

	ResultProcessor processor("moment.out");
	MomentEqSet eqSet;	RKSolver<double> solver(processor, eqSet);
	rk_params params;

	MomentEqSet::vec initialState(7);
	initialState = 0.0;
		params.final_time = 3e8;
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


	solver.solve(params, initialState);
}
