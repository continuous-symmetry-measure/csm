
#ifndef MOMENT_EQ_H
#define MOMENT_EQ_H

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "ChemicalNetwork.h"

class MomentEquations : public EquationSet<double>, public RKResultProcessor<double>, public ChemicalNetwork {
private:
	ofstream file;

public:


	virtual ~MomentEquations() {}
	MomentEquations(parsed_network pn) : ChemicalNetwork(pn) {}

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
				res[i] -= 2 * s.A *
				(state[chemicalTypes.size() + interactions.size() + interactionMat(i,i)] - state[i]);
			}

			// Check all interactions in which this is the input
			for (size_t j = 0; j < chemicalTypes.size(); ++j) {
				if (i != j && interactionMat(i,j) != -1) {
					int pos = interactionMat(i,j);
					const interaction &ii = interactions[pos];
					res[i] -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
							state[chemicalTypes.size() + pos];
				}
			}

			// Check all self interactions in which this is output (currently - only one...)
			for (size_t j = 0; j < indexedSelfOutputs[i].size(); ++j) {
				int pos = indexedSelfOutputs[i][j];
				const self_interaction &si = selfInteractions[pos];
				res[i] += chemicalTypes[si.input].A *
					(state[chemicalTypes.size() + interactions.size() + pos] - state[i]);
			}

			// Check all  interactions in which this is output (currently - only one...)
			for (size_t j = 0; j < indexedOutputs[i].size(); ++j) {
				int pos = indexedOutputs[i][j];
				interaction &ii = interactions[pos];
				res[i] += ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					(state[chemicalTypes.size() + pos]));
			}
		}

		// Go over all interactions
		for (size_t i = 0; i < interactions.size(); i++) {
			const interaction& ii = interactions[i];
			const species &s1 = chemicalTypes[ii.input1];
			const species &s2 = chemicalTypes[ii.input2];
			res[chemicalTypes.size() + i] =
				(s1.Flux * state[ii.input2] + s2.Flux * state[ii.input1] -
				 (s1.W + s2.W + s1.A + s2.A) * state[chemicalTypes.size() + i]);
		}

		for (size_t i = 0; i < selfInteractions.size(); i++) {
			const self_interaction &si = selfInteractions[i];
			const species &s = chemicalTypes[si.input];
			res[chemicalTypes.size() + interactions.size() + i] =
				s.Flux + 2 * s.Flux * state[si.input] + s.W * state[si.input] -
				2 * s.W	* state[chemicalTypes.size() + interactions.size() + i] -
				4 * s.A * (state[chemicalTypes.size() + interactions.size() + i] - state[si.input]);


			// Check all interactions in which this is the input
			for (size_t j = 0; j < chemicalTypes.size(); ++j) {
				if (si.input != j && interactionMat(si.input,j) != -1) {
					int pos = interactionMat(si.input, j);
					const interaction &ii = interactions[pos];
					res[chemicalTypes.size() + interactions.size() + i] -=
						(chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
							state[chemicalTypes.size() + pos];
				}
			}

			// Check all self interactions in which this is output (currently - only one...)
			for (size_t j = 0; j < indexedSelfOutputs[i].size(); ++j) {
				int pos = indexedSelfOutputs[i][j];
				const self_interaction &si2 = selfInteractions[pos];
				res[chemicalTypes.size() + interactions.size() + i] +=
					chemicalTypes[si2.input].A *
					(state[chemicalTypes.size() + interactions.size() + pos] - state[si2.input]);
			}

			// Check all  interactions in which this is output (currently - only one...)
			for (size_t j = 0; j < indexedOutputs[i].size(); ++j) {
				int pos = indexedOutputs[i][j];
				const interaction &ii = interactions[pos];
				res[chemicalTypes.size() + interactions.size() + i] +=
					((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
						(state[chemicalTypes.size() + pos]));
			}


		}

		return res;

	}

	vec createInitialConditions() {

		vec result(chemicalTypes.size() + interactions.size() + selfInteractions.size());
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
	virtual void stepPerformed(double time, double dt, const vec& state, const vec& prevState) {

		file << time << "\t" << dt << "\t";

		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			file << state[i] << "\t";
		}
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &ii = interactions[i];
			file << ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					state[i + chemicalTypes.size()]) << "\t";
		}

		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			file << (chemicalTypes[si.input].A *
					(state[i + chemicalTypes.size() + interactions.size()] - state[si.input])) << "\t";
		}

		file << endl;
	}

	/**
	 * Announce that the solving is complete
 	 */
	virtual void solutionComplete(const vec& finalState) {
		file << "Finished Solving" << endl;
		file.close();
	}
};

#endif
