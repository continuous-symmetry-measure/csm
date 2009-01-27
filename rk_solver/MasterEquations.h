#ifndef MASTER_EQ_H
#define MASTER_EQ_H

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "rk_solver.h"
#include "ChemicalNetwork.h"

using namespace std;

class MasterEquations : public EquationSet<double>, public RKResultProcessor<double>, public ChemicalNetwork{
private:
	size_t maxPos;
	vector<size_t> diffs;
	
	
	unsigned long stepNum;
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
	MasterEquations(parsed_network pn) : ChemicalNetwork(pn) {
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


			// Next, go over self interactions
			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];

				// Build the new index:
				size_t inputState = pos + 2 * diffs[si.input];
				bool use = true;
				for (size_t k = 0; k < si.outputs.size(); ++k) {
					if (si.outputs[k] < chemicalTypes.size()) {
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
					if (ii.outputs[k] < chemicalTypes.size()) {
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

			// Next, go over Dissociations
			for (size_t i = 0; i < dissociations.size(); i++) {
				const dissociation &di = dissociations[i];
				size_t i1 = indices[di.input];				

				// Build the new index:
				size_t inputState = pos + diffs[di.input];
				bool use = true;
				for (size_t k = 0; k < di.outputs.size(); ++k) {
					if (di.outputs[k] < chemicalTypes.size()) {
						if (indices[di.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[di.outputs[k]];
						}
					}
				}

				res[pos] += (di.D) *
					((((i1 == chemicalTypes[di.input].cutoff) || (!use)) ?
					0 : ((i1 + 1) * state[inputState])) - (i1 * state[pos]));
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

	/** 
	 * Compute Rn/Wn for iterative method. Rn - sum of changes to derivative entering from ext. sources. 
	 * wn - sum of coeff. of exiting state
	 * @param state
	 * @return Rn / Wn
	 */
	virtual vec comptueRnDivWn(const vec& state) { 
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
			double wn = 0.0;
			double rn = 0.0;

			// First - go over types
			for (size_t i = 0; i < chemicalTypes.size(); i++) {
				wn += (chemicalTypes[i].Flux + chemicalTypes[i].W * indices[i]);
				rn += chemicalTypes[i].Flux * ((indices[i] == 0) ? 0 : state[pos - diffs[i]])
					+ chemicalTypes[i].W * ((indices[i] == chemicalTypes[i].cutoff) ? 0 : (indices[i] + 1) * state[pos + diffs[i]]);
			}

			// Next, go over self interactions
			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];

				// Build the new index:
				size_t inputState = pos + 2 * diffs[si.input];
				bool use = true;
				for (size_t k = 0; k < si.outputs.size(); ++k) {
					if (si.outputs[k] < chemicalTypes.size()) {
						if (indices[si.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[si.outputs[k]];
						}
					}
				}

				wn += chemicalTypes[si.input].A * ii * (ii - 1);
				rn += chemicalTypes[si.input].A *
				    (((ii >= chemicalTypes[si.input].cutoff - 1) || (!use)) ?
				    		0 : ((ii + 2) * (ii + 1) * state[inputState]));
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
					if (ii.outputs[k] < chemicalTypes.size()) {
						if (indices[ii.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[ii.outputs[k]];
						}
					}
				}

				wn += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * i1 * i2;
				rn += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					(((i1 == chemicalTypes[ii.input1].cutoff) || (i2 == chemicalTypes[ii.input2].cutoff) || (!use)) ?
					 	0 : ((i1 + 1) * (i2 + 1) * state[inputState]));
			}

			// Next, go over Dissociations
			for (size_t i = 0; i < dissociations.size(); i++) {
				const dissociation &di = dissociations[i];
				size_t i1 = indices[di.input];				

				// Build the new index:
				size_t inputState = pos + diffs[di.input];

				bool use = true;
				for (size_t k = 0; k < di.outputs.size(); ++k) {
					if (di.outputs[k] < chemicalTypes.size()) {
						if (indices[di.outputs[k]] == 0) {
							use = false;
							break;
						} else {
							inputState -= diffs[di.outputs[k]];
						}
					}
				}
				
				wn += (di.D) * i1;
				rn += (di.D) *
					((i1 == chemicalTypes[di.input].cutoff) || (!use)) ?
					0 : ((i1 + 1) * state[inputState]);;
			}				
			
			res[pos] = rn / wn;

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
	virtual vec compute_derivative_derivative(const vec& state)  {
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

			// Go over self interactions
			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];
				res[pos] -= chemicalTypes[si.input].A * ii * (ii - 1);
			}

			// Next - Go over two-component interactions
			for (size_t i = 0; i < interactions.size(); i++) {
				const interaction &ii = interactions[i];
				size_t i1 = indices[ii.input1];
				size_t i2 = indices[ii.input2];
				res[pos] -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * i1 * i2;
			}

			// Go over dissociations
			for (size_t i = 0; i < dissociations.size(); i++) {
				const dissociation& di = dissociations[i];
				size_t i1 = indices[di.input];
				res[pos] -= di.D * i1;
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
    	 * Normalize the results (if relevant)
	 * 
	 * @param state The un-normalized state
	 * @return normalized result
	 */ 
	virtual vec normalize(const vec& state) { 
		return state.copy() / state.sum();
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
		
		for (size_t i = 0; i < dissociations.size(); i++) {
			dissociation &di = dissociations[i];
			for (size_t k = 0; k < di.outputs.size(); ++k) {
				file << "R_" << getOutputName(di.outputs[k]);
				if (k != di.outputs.size() - 1) {
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
	virtual void stepPerformed(double time, double dt, const vec& state, const vec& prevState, bool forcePrint = false) {

		stepNum++;
		if (stepNum % 40 != 1 && (!forcePrint)) return;

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

		for (size_t i = 0; i < dissociations.size(); i++) {
			dissociation &di = dissociations[i];
			file << di.D * (avgVec[di.input]) << "\t";
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
