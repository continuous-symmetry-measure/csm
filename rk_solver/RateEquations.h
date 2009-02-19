
#ifndef RATE_EQ_H
#define RATE_EQ_H

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "rk_solver.h"
#include "ChemicalNetwork.h"


class RateEquations : public EquationSet<double>, public RKResultProcessor<double>, public ChemicalNetwork {

private:
	ofstream file;

public:

	virtual ~RateEquations() {}
	RateEquations(parsed_network pn) : ChemicalNetwork(pn) {}

	/**
	 * Update any relevant parameters prior to the next step
	 *
	 * Default implementation is empty, assuming there are no varying parameters
	 *
	 * @param state The current (pre-step) state
	 */
	virtual void updateParameters(const vec& state) {}

	/**
 	 * Compute the error, given the previous state, the current state, and the used delta
  	 */
	virtual double computeError(const vec& state, const vec& prevState, double usedDelta) {
		return scalarMax(abs(state - prevState) / state) / usedDelta;
	}


	/**
	 * According to a current state, and the
	 * relevant equations, compute the time derivative.
	 *
	 * @param time The current time
	 * @param state The current state to compute derivative for
	 * @return The time derivative of the equations
	 */
	virtual vec compute_derivative(const double t, const vec& state) {
		vec res(state.size());
		res = 0;

		// First, go over the means of the populations
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			const species& s = chemicalTypes[i];
			res[i] += (s.Flux - s.W * state[i]);

			// Check if there is a self interaction for this:
			if (interactionMat(i,i) != -1) {
				res[i] -= 2 * s.A * (state[i] * state[i]);
			}

			// Check if there is a dissociation for this:
			if (dissociationMat(i,i) != -1) {
				res[i] -= dissociations[dissociationMat(i,i)].D * state[i];
			}

			// Check all interactions in which this is the input
			for (size_t j = 0; j < chemicalTypes.size(); ++j) {
				if (i != j && interactionMat(i,j) != -1) {
					int pos = interactionMat(i,j);
					const interaction &ii = interactions[pos];
					res[i] -= ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
							state[ii.input1] * state[ii.input2]);
				}
			}

			// Check all self interactions in which this is output
			for (size_t j = 0; j < indexedSelfOutputs[i].size(); ++j) {
				size_t pos = indexedSelfOutputs[i][j];
				const self_interaction &si = selfInteractions[pos];
				res[i] += chemicalTypes[si.input].A * state[si.input] * state[si.input];
			}

			// Check all  dissociations in which this is output 
			for (size_t j = 0; j < indexedDissociations[i].size(); ++j) {
				size_t pos = indexedDissociations[i][j];
				dissociation &di = dissociations[pos];
				res[i] += di.D * state[di.input];					
			}

			// Check all  interactions in which this is output 
			for (size_t j = 0; j < indexedOutputs[i].size(); ++j) {
				size_t pos = indexedOutputs[i][j];
				interaction &ii = interactions[pos];
				res[i] += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					(state[ii.input1] * state[ii.input2]);
			}
		}
		return res;
	}

	vec createInitialConditions() {
		vec result(chemicalTypes.size());
		result = 0;
		return result;
	}

	/**
	* Announce that the solving has started
	*/
	virtual void solvingStarted(const rk_params &params, const vec& initialState) {
		file.open("rate.out");
		file.width(10);
		file.precision(4);

		file << "Time\tdt\t";

		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			file << "<" << chemicalTypes[i].name << ">" << "\t\t";
		}

		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &ii = interactions[i];
			for (size_t k = 0; k < ii.outputs.size(); ++k) {
				file << "R_" << getOutputName(ii.outputs[k]);
				if (k != ii.outputs.size() - 1) {
					file << ",";
				} else {
					file << "\t";
				}
			}

		}

		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			for (size_t k = 0; k < si.outputs.size(); ++k) {
				file << "R_" << getOutputName(si.outputs[k]);
				if (k != si.outputs.size() - 1) {
					file << ",";
				} else {
					file << "\t";
				}
			}
		}
		file << endl;
	}

	/**
	 * A single complete RK step has been performed
 	 * @param time	The time
	 * @param dt	The chosen delta
	 * @param state The state after the step
	 */
	virtual void stepPerformed(double time, double dt, const vec& state, const vec& prevState, bool forceUpdate) {

		file << time << "\t" << dt << "\t";

		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			file << state[i] << "\t";
		}
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &ii = interactions[i];
			file << ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					state[ii.input1] * state[ii.input2]) << "\t";
		}

		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			file << chemicalTypes[si.input].A * state[si.input] * state[si.input] << "\t";
		}

		file << endl;
	}

	/**
	 * Announce that the solving is complete
 	 */
	virtual void solutionComplete(const vec& finalState) {
		file.close();
	}
};

#endif
