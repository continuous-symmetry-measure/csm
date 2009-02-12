#ifndef FULLPAIR_EQ
#define FULLPAIR_EQ

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "ChemicalNetwork.h"
#include "rk_solver.h"

// Here the state is a vector of matrices, each (cutoff1 + 1) x (cutoff2 + 1)
// sized according to the interaction types
class FullPairEquations : public ChemicalNetwork, public EquationSet<Matrix<double> >, public RKResultProcessor<Matrix<double> >{

protected:
	struct interaction_data {
		Vector<double> conditionalMean1;		// The mean values for the first input given the second value
		Vector<double> conditionalMean2;		// The mean values for the second input given the first value
		Vector<double> conditionalSecond1;		// The second moment of the first input given the second input
		Vector<double> conditionalSecond2;		// The second moment of the second input given the first
		double mean1;		// The total mean of the first input
		double mean2;		// The total mean of the second input
		double secondMoment1;		// The second moment
		double secondMoment2;		// The second moment
		double corr;			// <input1 * input2>		
	};

	Matrix<interaction_data> interactionData;			// The interactions interactionData

	ofstream file;

	unsigned long stepNum;

	Vector<double> prepareResultsVec() {
		Vector<double> results(5 * interactionData.nrows() * (interactionData.nrows() - 1) / 2);
		int pos = 0;
		for (size_t i = 0; i < interactionData.nrows(); i++) {
			for (size_t j = 0; j < i; j++) {
				interaction_data &inter = interactionData(i,j);
				results[pos++] = inter.mean1;
				results[pos++] = inter.mean2;
				results[pos++] = inter.secondMoment1;
				results[pos++] = inter.secondMoment2;
				results[pos++] = inter.corr;
			}
		}
		return results;
	}

	void extendInteractionInfo() {
		interactionData.resize(chemicalTypes.size(), chemicalTypes.size());
		for (size_t i = 0; i < interactionData.nrows(); i++) {
			for (size_t j = 0; j < i; j++) {
				interaction_data& id = interactionData(i,j);
				id.conditionalMean1.resize(chemicalTypes[j].cutoff + 1);
				id.conditionalMean2.resize(chemicalTypes[i].cutoff + 1);
				id.conditionalSecond1.resize(chemicalTypes[j].cutoff + 1);
				id.conditionalSecond2.resize(chemicalTypes[i].cutoff + 1);
				interactionData(j,i) = interactionData(i,j);
			}
		}
	}

public:

	virtual ~FullPairEquations() {}
	FullPairEquations(parsed_network pn) : ChemicalNetwork(pn) {
		extendInteractionInfo();
	}

	/**
	 * Update any relevant parameters prior to the next step
	 *
	 * Default implementation is empty, assuming there are no varying parameters
	 *
	 * @param state The current (pre-step) state
	 */
	virtual void updateParameters(const vec& state) {
		// Compute the mean values for each column and row of the interaction equations
		size_t pos = 0;
		for (size_t m = 0; m < interactionData.nrows(); m++) {
			for (size_t n = 0; n < m; n++) { 				
				interaction_data& id = interactionData(m,n);
				const Matrix<double> &probs = state[pos];
				pos++;
	
				// Compute means
				id.conditionalMean1   = 0.0;
				id.conditionalMean2   = 0.0;
				id.conditionalSecond1 = 0.0;
				id.conditionalSecond2 = 0.0;
				id.mean1 = 0.0;
				id.mean2 = 0.0;
				id.secondMoment1 = id.secondMoment2 = 0.0;
				id.corr    = 0.0;
				double rowSum = 0.0;
				double colSum = 0.0;
	
				for (size_t j = 0; j <= chemicalTypes[m].cutoff; j++) {
					rowSum = probs.row(j).sum();
					for (size_t k = 0; k <= chemicalTypes[n].cutoff; k++) {
						colSum = probs.column(k).sum();
						id.conditionalMean1[k] += (colSum > 0.0) ? (j * probs(j,k) / colSum) : 0.0;
						id.conditionalMean2[j] += (rowSum > 0.0) ? (k * probs(j,k) / rowSum) : 0.0;						
						id.conditionalSecond1[k] += (colSum > 0.0) ? (j * j * probs(j,k) / colSum) : 0.0;
						id.conditionalSecond2[j] += (rowSum > 0.0) ? (k * k * probs(j,k) / rowSum) : 0.0;							
						id.mean1 += j * probs(j,k);
						id.mean2 += k * probs(j,k);
						id.secondMoment1 += j * j * probs(j,k);
						id.secondMoment2 += k * k * probs(j,k);
						id.corr += k * j * probs(j,k);
					}
				}
			}
		}
	}

	/**
 	 * Compute the error, given the previous state, the current state, and the used delta
  	 */
	virtual double computeError(const vec& state, const vec& prevState, double usedDelta) {
		updateParameters(prevState);
		Vector<double> prev = prepareResultsVec();
		updateParameters(state);
		Vector<double> current = prepareResultsVec();
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
		updateParameters(state);
		vec res = state.copy();

		// First - go over the two-component interactions
		size_t iii = 0;
		for (size_t m = 0; m < interactionData.nrows(); m++) {
			for (size_t n = 0; n < m; n++) { 
				species& input1 = chemicalTypes[m];
				species& input2 = chemicalTypes[n];
				Matrix<double>& interRes = res[iii];
				const Matrix<double>& stateMat = state[iii];
				iii ++;
	

				for (size_t j = 0; j <= input1.cutoff; ++j) {
					for (size_t k = 0; k <= input2.cutoff; ++k) {

						// 1. Take care of the non interaction part and this pair interaction
						interRes(j, k) =
							// Non-interaction parts
							input1.Flux * ((j == 0 ? 0 : stateMat(j - 1, k))  - stateMat(j, k)) +
							input2.Flux * ((k == 0 ? 0 : stateMat(j, k - 1))  - stateMat(j, k)) +
							input1.W * (((j == input1.cutoff)? 0 : ((j + 1) * stateMat(j + 1, k))) - (j * stateMat(j, k))) +
							input2.W * (((k == input2.cutoff)? 0 : ((k + 1) * stateMat(j, k + 1))) - (k * stateMat(j, k)));
						
						// Go over interactions
						for (size_t interNum = 0; interNum < interactions.size(); interNum++) {
							interaction& inter = interactions[interNum];
							interaction_data& iid = interactionData(
								max(inter.input1, inter.input2), min(inter.input1, inter.input2));

							bool isFirstOutput = false;
							bool isSecondOutput = false;
							for (size_t p = 0; p < inter.outputs.size(); p++) {
								if (m == inter.outputs[p]) isFirstOutput = true;
								if (n == inter.outputs[p])  isSecondOutput = true;
							}

							// if this interaction is between the two species:
							if ((inter.input1 == m && inter.input2 == n) || (inter.input1 == n && inter.input2 == m)) {
								interRes(j,k) += 						
									(input1.A + input2.A) * (((j == input1.cutoff) || (k == input2.cutoff) ? 0 :
									((j+1)*(k+1)*stateMat(j+1, k+1))) - (j * k * stateMat(j,k)));
							} else if (inter.input1 == m || inter.input2 == m) {
								// The m species is between the reactants								
								Vector<double> &means = (m == inter.input1) ? iid.conditionalMean2 : iid.conditionalMean1;
								if (isSecondOutput) {
									interRes(j, k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) *
										((((j == input1.cutoff) || k == 0)? 0 : 
										(j + 1) * means[j+1] * stateMat(j+1, k - 1)) -
										(j * means[j] * stateMat(j,k)));
								} else {
									interRes(j, k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) *
										(((j == input1.cutoff) ? 0 : 
										(j + 1) * means[j+1] * stateMat(j+1, k)) -
										(j * means[j] * stateMat(j,k)));
								}
							} else if (inter.input1 == n || inter.input2 == n) {
								// The n species is between the reactants								
								Vector<double> &means = (n == inter.input1) ? iid.conditionalMean2 : iid.conditionalMean1;
								if (isFirstOutput) {
									interRes(j, k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) *
										((((k == input1.cutoff) || j == 0)? 0 : 
										(k + 1) * means[k+1] * stateMat(j - 1, k + 1)) -
										(k * means[k] * stateMat(j,k)));
								} else {
									interRes(j, k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) *
										(((k == input1.cutoff) ? 0 : 
										(k + 1) * means[k+1] * stateMat(j, k + 1)) -
										(k * means[k] * stateMat(j,k)));
								}
							} else {
								// Neither is one of the reactants
								if (isFirstOutput && isSecondOutput) {
									interRes(j,k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) * 
										iid.corr * (((j == 0 || k == 0) ? 0 : stateMat(j - 1, k - 1)) - stateMat(j,k));
								} else if (isFirstOutput) {
									interRes(j,k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) * 
										iid.corr * (((j == 0) ? 0 : stateMat(j - 1, k)) - stateMat(j,k));
								} else if (isSecondOutput) {
									interRes(j,k) += (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A) * 
										iid.corr * (((k == 0) ? 0 : stateMat(j, k - 1)) - stateMat(j,k));
								}
							}
						}

						// Go over self interactions
						for (size_t selfInterNum = 0; selfInterNum < selfInteractions.size(); selfInterNum++) {
							self_interaction selfInter = selfInteractions[selfInterNum];
							
							bool isFirstOutput = false;
							bool isSecondOutput = false;
							for (size_t p = 0; p < selfInter.outputs.size(); p++) {
								if (selfInter.outputs[p] == m) isFirstOutput = true;
								if (selfInter.outputs[p] == n) isSecondOutput = true;
							}

							// Check if the m species is the input
							if (selfInter.input == m) {
								if (isSecondOutput) {
									interRes(j,k) += chemicalTypes[selfInter.input].A *
										(((((j + 2) <= input1.cutoff) && (k > 0)) ? (j+2)*(j+1) * stateMat(j + 2, k - 1) : 0 ) - 
											(j * (j - 1) * stateMat(j,k)));
								} else {
									interRes(j,k) += chemicalTypes[selfInter.input].A *
										((((j + 2) <= input1.cutoff) ? (j+2)*(j+1) * stateMat(j + 2, k) : 0 ) - 
											(j * (j - 1) * stateMat(j,k)));
								}
							} else if (selfInter.input == n) {
								// Check if the n species is the input
								if (isFirstOutput) {
									interRes(j,k) += chemicalTypes[selfInter.input].A *
										(((((k + 2) <= input2.cutoff) && (j > 0)) ? (k+2)*(k+1) * stateMat(j - 1, k + 2) : 0 ) - 
										(k * (k - 1) * stateMat(j,k)));
								} else {
									interRes(j,k) += chemicalTypes[selfInter.input].A *
										((((k + 2) <= input2.cutoff) ? (k+2)*(k+1) * stateMat(j, k + 2) : 0 ) - 
										(k * (k - 1) * stateMat(j,k)));
								}
							} else {
								// They only appear as outputs
								// if both are outputs, we shall pick the second moment - first moment from some interaction, 
								// unconditioned
								if (isFirstOutput && isSecondOutput) {
									interaction_data& iid = interactionData(max(selfInter.input, m), min(selfInter.input,m));
									double param = 0;
									if (selfInter.input > m) {
										param = iid.secondMoment1 - iid.mean1;
									} else {
										param = iid.secondMoment2 - iid.mean2;
									}
									interRes(j,k) += chemicalTypes[selfInter.input].A * param * 
										(((j == 0 || k == 0) ? 0 : stateMat(j - 1, k - 1)) - stateMat(j,k));										
								} else if (isFirstOutput) {
									interaction_data& iid = interactionData(max(selfInter.input, m), min(selfInter.input,m));
									double param1 = 0;
									double param2 = 0;
									if (selfInter.input > m) {
										if (j != 0) {
											param1 = iid.conditionalSecond1[j - 1] - iid.conditionalMean1[j-1]; 
										}
										param2 = iid.conditionalSecond1[j] - iid.conditionalMean1[j];
									} else {
										if (j != 0) {
											param1 = iid.conditionalSecond2[j - 1] - iid.conditionalMean2[j-1]; 
										}
										param2 = iid.conditionalSecond2[j] - iid.conditionalMean2[j];
									}
									interRes(j,k) += chemicalTypes[selfInter.input].A * 
										(((j == 0) ? 0 : param1 * stateMat(j - 1, k)) - param2 * stateMat(j,k));
								} else if (isSecondOutput) {
									interaction_data& iid = interactionData(max(selfInter.input, n), min(selfInter.input,n));
									double param1 = 0;
									double param2 = 0;
									if (selfInter.input > n) {
										if (k != 0) {
											param1 = iid.conditionalSecond1[k - 1] - iid.conditionalMean1[k-1]; 
										}
										param2 = iid.conditionalSecond1[k] - iid.conditionalMean1[k];
									} else {
										if (k != 0) {
											param1 = iid.conditionalSecond2[k - 1] - iid.conditionalMean2[k-1]; 
										}
										param2 = iid.conditionalSecond2[k] - iid.conditionalMean2[k];
									}
									interRes(j,k) += chemicalTypes[selfInter.input].A * 
										(((k == 0) ? 0 : param1 * stateMat(j, k - 1)) - param2 * stateMat(j,k));
								}
							}
						}
						// Go over self interactions
						for (size_t disNum = 0; disNum < dissociations.size(); disNum++) {
							dissociation dis = dissociations[disNum];

							bool isFirstOutput = false;
							bool isSecondOutput = false;
							for (size_t p = 0; p < dis.outputs.size(); p++) {
								if (dis.outputs[p] == m) isFirstOutput = true;
								if (dis.outputs[p] == n) isSecondOutput = true;
							}

							// Check if the m species is the input
							if (dis.input == m) {
								if (isSecondOutput) {
									interRes(j,k) += dis.D * 
										(((((j + 1) <= input1.cutoff) && (k > 0)) ? (j+1) * stateMat(j + 1, k - 1) : 0 ) - 
										(j * stateMat(j,k)));
								} else {
									interRes(j,k) += dis.D * 
										((((j + 1) <= input1.cutoff) ? (j+1) * stateMat(j + 1, k) : 0 ) - 
										(j * stateMat(j,k)));
								}
							} else if (dis.input == n) {
								// Check if the n species is the input
								if (isFirstOutput) {
									interRes(j,k) += dis.D * 
										(((((k + 1) <= input2.cutoff) && (j > 0)) ? (k+1) * stateMat(j - 1, k + 1) : 0 ) - 
										(k * stateMat(j,k)));
								} else {
									interRes(j,k) += dis.D * 
										((((k + 1) <= input2.cutoff) ? (k+1) * stateMat(j, k + 1) : 0 ) - 
										(k * stateMat(j,k)));
								}
							} else {
								// They only appear as outputs
								// if both are outputs, we shall pick the second moment - first moment from some interaction, 
								// unconditioned
								if (isFirstOutput && isSecondOutput) {
									interaction_data& iid = interactionData(max(dis.input, m), min(dis.input,m));
									double param = 0;
									if (dis.input > m) {
										param = iid.mean1;
									} else {
										param = iid.mean2;
									}
									interRes(j,k) += dis.D * param * 
										(((j == 0 || k == 0) ? 0 : stateMat(j - 1, k - 1)) - stateMat(j,k));										
								} else if (isFirstOutput) {
									interaction_data& iid = interactionData(max(dis.input, m), min(dis.input,m));
									double param1 = 0;
									double param2 = 0;
									if (dis.input > m) {
										if (j != 0) {
											param1 = iid.conditionalMean1[j-1]; 
										}
										param2 = iid.conditionalMean1[j];
									} else {
										if (j != 0) {
											param1 = iid.conditionalMean2[j-1]; 
										}
										param2 = iid.conditionalMean2[j];
									}
									interRes(j,k) += dis.D * 
										(((j == 0) ? 0 : param1 * stateMat(j - 1, k)) - param2 * stateMat(j,k));
								} else if (isSecondOutput) {
									interaction_data& iid = interactionData(max(dis.input, n), min(dis.input,n));
									double param1 = 0;
									double param2 = 0;
									if (dis.input > n) {
										if (k != 0) {
											param1 = iid.conditionalMean1[k-1]; 
										}
										param2 = iid.conditionalMean1[k];
									} else {
										if (k != 0) {
											param1 = iid.conditionalMean2[k-1]; 
										}
										param2 = iid.conditionalMean2[k];
									}
									interRes(j,k) += dis.D * 
										(((k == 0) ? 0 : param1 * stateMat(j, k - 1)) - param2 * stateMat(j,k));
								}
							}
						}
					}
				}
			}
		}	
		return res;
	}

	/** 
	* For newton's method of finding steady state - 
	* compute the partial derivative of the equation for the time derivative
	* according the the variable it belongs to
	* 
	* @param state The current state to compute derivative for
	* @return the derivative's partial derivative
	*/

	vec createInitialConditions() {
		vec result(interactionData.nrows() * (interactionData.nrows() - 1) / 2);
		// Start with empty state;
		size_t iii = 0;
		for (size_t i = 0; i < interactionData.nrows(); i++) {
			for (size_t j = 0; j < i; j++) {
				result[iii].resize(
					chemicalTypes[i].cutoff + 1,
					chemicalTypes[j].cutoff + 1);				
				result[iii] = 0.0;
				result[iii](0,0) = 1.0;
				iii ++;
			}
		}
		return result;
	}

	/**
	* Announce that the solving has started
	*/
	virtual void solvingStarted(const rk_params &params, const vec& initialState) {
		file.open("pairs.out");
		file.width(10);
		file.precision(4);

		file << "Time\tdt\t";
		for (size_t i = 0; i < interactions.size(); ++i) {
			interaction &ii = interactions[i];
			file << "<" << chemicalTypes[ii.input1].name << ">" << "\t";
			file << "<" << chemicalTypes[ii.input2].name << ">" << "\t";
			for (size_t k = 0; k < ii.outputs.size(); k++) {
				file << "R_" << getOutputName(ii.outputs[k]);
				if (k != ii.outputs.size() - 1) {
					 file << ",";
				} else {
					file << "\t";
				}
			}

			if (interactionMat(ii.input1,ii.input1) != -1) {
				self_interaction &si = selfInteractions[interactionMat(ii.input1,ii.input1)];
				for (size_t k = 0; k < si.outputs.size(); k++) {
					file << "R_" << getOutputName(si.outputs[k]);
					if (k != si.outputs.size() - 1) {
						 file << ",";
					} else {
						file << "\t";
					}
				}
			}
			if (interactionMat(ii.input2,ii.input2) != -1) {
				self_interaction &si = selfInteractions[interactionMat(ii.input2,ii.input2)];
				for (size_t k = 0; k < si.outputs.size(); k++) {
					file << "R_" << getOutputName(si.outputs[k]);
					if (k != si.outputs.size() - 1) {
						 file << ",";
					} else {
						file << "\t";
					}
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
		updateParameters(state);
		const double MAX_DEV = 0.0001;
		if (!forcePrint) {
			for (size_t i = 0; i < interactions.size(); i++) {
				const Matrix<double> &probs = state[i];

				// Sum over all just to test...
				double sum = 0;
				for (size_t j = 0; j < probs.nrows(); j++) {
					sum += probs.row(j).sum();
				}
				if (fabs(sum - 1.0) > MAX_DEV) {
					cout << "Probability does not sum up to 1: " << sum << " - exiting..." << endl;	
					exit(1);
				}
			}
			stepNum++;
			if (stepNum % 50 != 1) return;
		}

		file << time << "\t" << dt << "\t";
		for (size_t i = 0; i < interactions.size(); ++i) {
			const interaction& inter = interactions[i];
			const interaction_data &id = interactionData(max(inter.input1, inter.input2), min(inter.input1,inter.input2));
			file << id.mean1 << "\t";
			file << id.mean2 << "\t";
			file << (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A)*(id.corr) << "\t";
			if (interactionMat(inter.input1,inter.input1) != -1) {
				file << chemicalTypes[inter.input1].A * (id.secondMoment1 - id.mean1) << "\t";
			}
			if (interactionMat(inter.input2,inter.input2) != -1) {
				file << chemicalTypes[inter.input2].A * (id.secondMoment2 - id.mean2) << "\t";
			}
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
