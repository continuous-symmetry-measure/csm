
#ifndef MASTER_EQ_H
#define MASTER_EQ_H

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "ChemicalNetwork.h"

using namespace std;

class MasterEquations : public EquationSet<double>, public RKResultProcessor<double>, public ChemicalNetwork{
private:
	size_t maxPos;
	vector<size_t> diffs;

	ofstream file;

	// Currently, computes only avarages
	vec prepareResultsVec(const vec& state) {
		vec results(chemicalTypes.size());
		results = 0;

		// Generic code for enumeration
		vector<size_t> indices(chemicalTypes.size());
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			indices[i] = 0;
		}
		size_t pos = 0;

		while (pos < maxPos) {
			// Operate on current state
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				results[i] += indices[i] * state[pos];
			}

			// Advance...
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				indices[i] += 1;
				if (indices[i] == chemicalTypes[i].cutoff + 1) {
					indices[i] = 0;
				} else {
					break;
				}
			}
			pos++;
		}
		return results;
	}

	// Currently, computes only avarages
	vec prepareSecondMomentVec(const vec& state) {
		vec results(chemicalTypes.size());
		results = 0;

		// Generic code for enumeration
		vector<size_t> indices(chemicalTypes.size());
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			indices[i] = 0;
		}
		size_t pos = 0;

		while (pos < maxPos) {
			// Operate on current state
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				results[i] += indices[i] * indices[i] * state[pos];
			}

			// Advance...
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				indices[i] += 1;
				if (indices[i] == chemicalTypes[i].cutoff + 1) {
					indices[i] = 0;
				} else {
					break;
				}
			}
			pos++;
		}
		return results;
	}

	// Prepare Correlation Vector
	vec prepareCorrVec(const vec& state) {
		vec results(interactions.size());
		results = 0;

		// Generic code for enumeration
		vector<size_t> indices(chemicalTypes.size());
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			indices[i] = 0;
		}
		size_t pos = 0;

		while (pos < maxPos) {
			// Operate on current state
			for (size_t i = 0; i < interactions.size(); ++i) {
				results[i] += indices[interactions[i].input1] * indices[interactions[i].input2] * state[pos];
			}

			// Advance...
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				indices[i] += 1;
				if (indices[i] == chemicalTypes[i].cutoff + 1) {
					indices[i] = 0;
				} else {
					break;
				}
			}
			pos++;
		}
		return results;
	}


public:

	virtual ~MasterEquations() {}
	MasterEquations(parsed_network pn) : ChemicalNetwork(pn) {}

	virtual void extendInteractionInfo() {
		maxPos = 1;
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			maxPos *= (chemicalTypes[i].cutoff + 1);
			if (i == 0) {
				diffs.push_back(1);
			} else {
				diffs.push_back(diffs[i - 1] * (chemicalTypes[i - 1].cutoff + 1));
			}
		}
	}

	/**
 	 * Compute the error, given the previous state, the current state, and the used delta
  	 */
	virtual double computeError(const vec& state, const vec& prevState, double usedDelta) {
		vec prev = prepareResultsVec(prevState);
		vec current = prepareResultsVec(state);
		return scalarMax(abs(current - prev) / prev) / usedDelta;
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
		vec res(maxPos);
		res = 0;

		// Generic code for enumeration
		vector<size_t> indices(chemicalTypes.size());
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			indices[i] = 0;
		}
		size_t pos = 0;

		while (pos < maxPos) {
			// Operate on current state

			// First - go over types
			for (size_t i = 0; i < chemicalTypes.size(); i++) {
				res[pos] +=
					// Flux Term
					chemicalTypes[i].Flux * (((indices[i] == 0) ? 0 : state[pos - diffs[i]]) - state[pos]) +
					// W Term
					chemicalTypes[i].W *
					(((indices[i] == chemicalTypes[i].cutoff) ? 0 : (indices[i] + 1) * state[pos + diffs[i]]) -
					(indices[i] * state[pos]));
			}


			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];

				// Build the new index:
				size_t inputState = pos + 2 * diffs[si.input];
				bool use = true;
				for (size_t k = 0; k < si.outputs.size(); ++k) {
					if (si.outputs[k] >= chemicalTypes.size()) {
						if (indices[si.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[si.outputs[k]];
						}
					}
				}

				res[pos] += chemicalTypes[si.input].A *
				    ((((ii >= chemicalTypes[si.input].cutoff - 1) || (!use)) ?
				    		0 : ((ii + 2) * (ii + 1) * state[inputState])) -
				    		(ii * (ii - 1) * state[pos]));
			}

			// Next - Go over two-component interactions
			for (size_t i = 0; i < interactions.size(); i++) {
				const interaction &ii = interactions[i];
				size_t i1 = indices[ii.input1];
				size_t i2 = indices[ii.input2];

				// Build the new index:
				size_t inputState = pos + diffs[ii.input1] + diffs[ii.input2];
				bool use = true;
				for (size_t k = 0; k < ii.outputs.size(); ++k) {
					if (ii.outputs[k] >= chemicalTypes.size()) {
						if (indices[ii.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[ii.outputs[k]];
						}
					}
				}

				res[pos] += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					((((i1 == chemicalTypes[ii.input1].cutoff) || (i2 == chemicalTypes[ii.input2].cutoff) || (!use)) ?
					 	0 : ((i1 + 1) * (i2 + 1) * state[inputState])) -
					(i1 * i2 * state[pos]));
			}


			// Advance...
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				indices[i] += 1;
				if (indices[i] == chemicalTypes[i].cutoff + 1) {
					indices[i] = 0;
				} else {
					break;
				}
			}
			pos++;
		}
		return res;

	}

	// Compute the partial derivative of each derivative equation according to the probability variable
	// used for 1D newton method
	vec compute_derivative_derivative(const vec& state)  {
		vec res(maxPos);
		res = 0;

		// Generic code for enumeration
		vector<size_t> indices(chemicalTypes.size());
		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			indices[i] = 0;
		}
		size_t pos = 0;

		while (pos < maxPos) {
			// Operate on current state

			// First - go over types
			for (size_t i = 0; i < chemicalTypes.size(); i++) {
				res[pos] -=
					// Flux Term
					chemicalTypes[i].Flux +
					// W Term
					chemicalTypes[i].W * indices[i];
			}


			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];
				res -= chemicalTypes[si.input].A * ii * (ii - 1);
			}

			// Next - Go over two-component interactions
			for (size_t i = 0; i < interactions.size(); i++) {
				const interaction &ii = interactions[i];
				size_t i1 = indices[ii.input1];
				size_t i2 = indices[ii.input2];
				res -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * (i1 * i2);
			}


			// Advance...
			for (size_t i = 0; i < chemicalTypes.size(); ++i) {
				indices[i] += 1;
				if (indices[i] == chemicalTypes[i].cutoff + 1) {
					indices[i] = 0;
				} else {
					break;
				}
			}
			pos++;
		}
		return res;
	}

	vec createInitialConditions() {

		vec result(maxPos);
		result = 0;
		result[0] = 1.0;

		return result;
	}

	/**
 	 * Perform newton's method to find Steady State Solutions
	 */
	void performNewtonMethod(double max_err) {
		vec state = createInitialConditions();
		double err = max_err + 1;
		while (err > max_err) {
			vec newState = state - (compute_derivative(0.0,state) / compute_derivative_derivative(state));
			err = abs((newState - state) / state).max();
			state = newState;
			cout << "Error: " << err << endl;
		}

		file.open("master_newton.out");
		vec avgVec = prepareResultsVec(state);
		vec momentVec = prepareSecondMomentVec(state);
		vec corrVec = prepareCorrVec(state);

		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			file << "Mean of " << chemicalTypes[i].name << " is " << avgVec[i] << endl;
		}
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			file << "Production Rate of ";
			for (size_t k = 0; k < si.outputs.size(); ++k) {
				file << "R_" << getOutputName(si.outputs[k]);
				if (k != si.outputs.size() - 1) {
					file << ",";
				}
			}
			file << " is " << (chemicalTypes[si.input].A * (momentVec[si.input] - avgVec[si.input])) << endl;
		}

		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &ii = interactions[i];
			file << "Production Rate of ";
			for (size_t k = 0; k < ii.outputs.size(); ++k) {
				file << "R_" << getOutputName(ii.outputs[k]);
				if (k != ii.outputs.size() - 1) {
					file << ",";
				}
			}
			file << " is " << ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input1].A) * corrVec[i]) << endl;
		}
		file.close();

	}


	/**
	* Announce that the solving has started
	*/
	virtual void solvingStarted(const rk_params &params, const vec& initialState) {
		file.open("master.out");
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
		vec avgVec = prepareResultsVec(state);
		vec momentVec = prepareSecondMomentVec(state);
		vec corrVec = prepareCorrVec(state);

		file << time << "\t" << dt << "\t";

		for (size_t i = 0; i < chemicalTypes.size(); i++) {
			file << avgVec[i] << "\t";
		}
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &ii = interactions[i];
			file << ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input1].A) * corrVec[i])
				<< "\t";
		}

		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			file << (chemicalTypes[si.input].A * (momentVec[si.input] - avgVec[si.input])) << "\t";
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
