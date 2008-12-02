/**
 * A simple solver for a system of linear equations using
 * rk4-5
 */
#ifndef _RK_SOLVER_H_
#define _RK_SOLVER_H_

#include "equation_set.h"

struct rk_params {
	double initialDelta;	// The initial delta to try
	double max_error;	// The maximal allowed error;
	double initial_time;	// The initial time
	double final_time;	// The final time	
	double limit;		// The maximum change in the derivative which constitutes steady-state
};

template <typename T>
class RKResultProcessor {
public:

	typedef typename EquationSet<T>::vec res_vec;

	/**
	 * Announce that the solving has started
	 */
	virtual void solvingStarted(const rk_params &params, const res_vec& initalState) { }

	/** 
	 * A single complete RK step has been performed
	 * @param time	The time
	 * @param dt	The chosen delta
	 * @param state The state after the step
	 */
	virtual void stepPerformed(double time, double dt, const res_vec& state, const res_vec& prevState) { }
	
	/** 
	 * Announce that the solving is complete
	 */
	virtual void solutionComplete(const res_vec& finalState) { }

	virtual ~RKResultProcessor() {}
};
template <typename T>
class RKSolver {
public:
	typedef typename EquationSet<T>::vec vec;
private:
	rk_params params;
	RKResultProcessor<T>& processor;
	EquationSet<T>& eqSet;
	vec state; 
	double time;
	double delta;
	double usedDelta;

	// Perform a single RK step
	vec rk_step(double time, double delta, const vec& data) {
		vec k1 = eqSet.compute_derivative(time, data);		
		vec k2 = eqSet.compute_derivative(time + 
			delta / 2.0, data.copy() + (k1 * (delta * 0.5)));
		vec k3 = eqSet.compute_derivative(time + 
			delta / 2.0, data.copy() + (k2 * (delta * 0.5)));		
		vec k4 = eqSet.compute_derivative(time + delta, data.copy() + (k3 * delta));
		vec step = (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * ((1.0/6.0) * delta);
		return data.copy() + step;
	}

	// Compute the (non-relative) error
	double computeError(const vec& v1, const vec& v2) {
		return scalarMax(abs(v1 - v2));
	}

	vec performStep() {
		double error = params.max_error + 1;
		bool firstTime = true;
		vec oneStep(state.size()), twoStep(state.size());
		while (error > params.max_error) {
			if (firstTime) {
				firstTime = false;
			} else {
				delta = 0.9 * delta * 
					pow(params.max_error / error, 0.25);
			}
			// Advance one step
			oneStep = rk_step(time, delta, state);

			// Advance two steps
			twoStep = rk_step(time, delta / 2.0, state);
			twoStep = rk_step(time + delta / 2.0, delta / 2.0, 
				twoStep.copy());

			// Compute error between them
			error = computeError(oneStep, twoStep);
		}	
		time = time + delta;
		usedDelta = delta;
		
		double tempDelta;
		if (error > 0) {
			tempDelta = 0.9 * delta * pow((params.max_error / error), 0.2);
		} else {
			tempDelta = 3.0 * delta;
		}
		if (tempDelta > (3 * delta) || (tempDelta / delta) < 1.01) {
			delta *= 3;
		} else {
			delta = tempDelta;
		}

		return oneStep;
	}	

public: 

	RKSolver(RKResultProcessor<T>& resProcessor, EquationSet<T>& equations) : processor(resProcessor), eqSet(equations) {}

	void solve(const rk_params& paramSet, vec& initialState) {

		// Copy Params
		params = paramSet;

		// Initialize
		state.reshape(initialState.size());
		state = initialState;
		time = params.initial_time;
		delta = params.initialDelta;
		double error = params.limit + 1;

		processor.solvingStarted(params, initialState);				

		while (time < params.final_time && error > params.limit) {
			vec result = performStep();
			error = eqSet.computeError(result, state, usedDelta);
			processor.stepPerformed(time, delta, result, state);
			state = result;
			std::cout << "Error: " << error << std::endl;
		}
		
		processor.solutionComplete(state);
	}
};
#endif
