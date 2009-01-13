#ifndef PAIR_EQ
#define PAIR_EQ

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>
#include "ChemicalNetwork.h"
#include "rk_solver.h"

// Here the state is a vector of matrices, each (cutoff1 + 1) x (cutoff2 + 1)
// sized according to the interaction types
class PairEquations : public ChemicalNetwork, public EquationSet<Matrix<double> >, public RKResultProcessor<Matrix<double> >{

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

	vector<interaction_data> interactionData;			// The interactions interactionData

	ofstream file;

	unsigned long stepNum;

	Vector<double> prepareResultsVec() {
		Vector<double> results(5 * interactions.size());
		int pos = 0;
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction_data &inter = interactionData[i];
			results[pos++] = inter.mean1;
			results[pos++] = inter.mean2;
			results[pos++] = inter.secondMoment1;
			results[pos++] = inter.secondMoment2;
			results[pos++] = inter.corr;
		}
		return results;
	}

	void extendInteractionInfo() {
		for (size_t i = 0; i < interactions.size(); i++) {
			const interaction &ii = interactions[i];
			interaction_data id;
			id.conditionalMean1.resize(chemicalTypes[ii.input2].cutoff + 1);
			id.conditionalMean2.resize(chemicalTypes[ii.input1].cutoff + 1);
			id.conditionalSecond1.resize(chemicalTypes[ii.input2].cutoff + 1);
			id.conditionalSecond2.resize(chemicalTypes[ii.input1].cutoff + 1);
			interactionData.push_back(id);
		}
	}

public:

	virtual ~PairEquations() {}
	PairEquations(parsed_network pn) : ChemicalNetwork(pn) {
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
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction& inter = interactions[i];
			interaction_data& id = interactionData[i];
			const Matrix<double> &probs = state[i];

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

			for (size_t j = 0; j <= chemicalTypes[inter.input1].cutoff; j++) {
				rowSum = probs.row(j).sum();
				for (size_t k = 0; k <= chemicalTypes[inter.input2].cutoff; k++) {
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
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction& inter = interactions[i];
			species& input1 = chemicalTypes[inter.input1];
			species& input2 = chemicalTypes[inter.input2];
			Matrix<double>& interRes = res[i];

			const Matrix<double>& stateMat = state[i];


			for (size_t j = 0; j <= input1.cutoff; ++j) {
				for (size_t k = 0; k <= input2.cutoff; ++k) {

					// 1. Take care of the non interaction part and this pair interaction
					interRes(j, k) =
						// Non-interaction parts
						input1.Flux * ((j == 0 ? 0 : stateMat(j - 1, k))  - stateMat(j, k)) +
						input2.Flux * ((k == 0 ? 0 : stateMat(j, k - 1))  - stateMat(j, k)) +
						input1.W * (((j == input1.cutoff)? 0 : ((j + 1) * stateMat(j + 1, k))) - (j * stateMat(j, k))) +
						input2.W * (((k == input2.cutoff)? 0 : ((k + 1) * stateMat(j, k + 1))) - (k * stateMat(j, k))) +
						// current pair interaction
						(input1.A + input2.A) * (((j == input1.cutoff) || (k == input2.cutoff) ? 0 :
							((j+1)*(k+1)*stateMat(j+1, k+1))) - (j * k * stateMat(j,k)));

					// 2. Go over self interactions:
					// A. Self interactions where these are the inputs
					if (interactionMat(inter.input1, inter.input1) != -1) {
						// Check if the other is one of the outputs
						if (selfInteractionMat(inter.input1, inter.input2) != -1) {
							interRes(j,k) += input1.A *
								(((((j + 2) <= input1.cutoff) && k > 0) ? 
									((j+1)*(j+2)*stateMat(j+2, k - 1)) : 0) - (j * (j-1)*stateMat(j,k)));
						} else {
							interRes(j,k) += input1.A *
								((((j + 2) <= input1.cutoff) ? 
									((j+1)*(j+2)*stateMat(j+2, k)) : 0) - (j * (j-1)*stateMat(j,k)));
						}
					}
					if (interactionMat(inter.input2, inter.input2) != -1) {
						// Check if the other is one of the outputs					
						if (selfInteractionMat(inter.input2, inter.input1) != -1) {
							interRes(j,k) += input2.A *
								(((((k + 2) <= input2.cutoff) && j > 0)? 
									((k+1)*(k+2)*stateMat(j - 1, k+2)) : 0) - (k * (k-1)*stateMat(j,k)));
						} else {
							interRes(j,k) += input2.A *
								((((k + 2) <= input2.cutoff) ? 
									((k+1)*(k+2)*stateMat(j, k+2)) : 0) - (k * (k-1)*stateMat(j,k)));
						}
					}

					// B. Self interactions with either of these is the output
					// There are several cases - both are the output of the self interaction, one is, or none are
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l == inter.input1 || l == inter.input2) { 
							// One of the interacting types is the input, continue since we've taken care of this
							continue;
						} 
						int first = selfInteractionMat(l, inter.input1);
						int second = selfInteractionMat(l, inter.input2);
						int firstIn = interactionMat(l, inter.input1);
						int secondIn = interactionMat(l, inter.input2);

						double param1 = 0;
						double param2 = 0;
						if (first != -1 || second != -1) {													
							// If the reactants produces only one and does not interact with it, 
							// or produces both and either does not interact with none
							// interactes with both (special case)
							if (((first != -1 && second != -1) && 
									((firstIn != -1 && secondIn != -1) || (firstIn == -1 && secondIn == -1))) ||
									(first != -1 && firstIn == -1) || (second != -1 && secondIn == -1)) {	
								interaction &ii = interactions[selfInteractions[first].locations[0]];
								interaction_data &iid = interactionData[selfInteractions[first].locations[0]];
								// Check if it is the first or second input
								param1 = param2 = (l == ii.input1) ?
									iid.secondMoment1 - iid.mean1 : 
									iid.secondMoment2 - iid.mean2; 							
								} else if (first != -1 && firstIn != -1) {
									interaction &ii = interactions[interactionMat(l, inter.input1)];
									interaction_data &iid = interactionData[interactionMat(l, inter.input1)];
									// Check if it is the first or second input								
									if (j != 0) {
										param1 = (l == ii.input1) ?
										iid.conditionalSecond1[j - 1] - iid.conditionalMean1[j - 1] : 
										iid.conditionalSecond2[j - 1] - iid.conditionalMean2[j - 1]; 								
									}
									param2 = (l == ii.input1) ?
										iid.conditionalSecond1[j] - iid.conditionalMean1[j] : 
										iid.conditionalSecond2[j] - iid.conditionalMean2[j]; 							
								} else if (second != -1 && secondIn != -1) {
									interaction &ii = interactions[interactionMat(l, inter.input2)];
									interaction_data &iid = interactionData[interactionMat(l, inter.input2)];
									// Check if it is the first or second input								
									if (k != 0) {
										param1 = (l == ii.input1) ?
											iid.conditionalSecond1[k - 1] - iid.conditionalMean1[k - 1] : 
											iid.conditionalSecond2[k - 1] - iid.conditionalMean2[k - 1]; 								
									}
									param2 = (l == ii.input1) ?
										iid.conditionalSecond1[k] - iid.conditionalMean1[k] : 
										iid.conditionalSecond2[k] - iid.conditionalMean2[k]; 							
								} else {
									vmerror("BUGBUGBUG");
								}
							}
							if (first != -1 && second != -1) { 
								interRes(j,k) += chemicalTypes[l].A * 
									(param1 * (j == 0 || k == 0 ? 0 : stateMat(j - 1, k - 1)) - param2 * stateMat(j, k));
							} else if (first != -1) {
								interRes(j,k) += chemicalTypes[l].A * 
									(param1 * (j == 0 ? 0 : stateMat(j - 1, k)) - param2 * stateMat(j, k));
							} else if (second != -1) { 
								interRes(j,k) += chemicalTypes[l].A *
								(param1 * (k == 0 ? 0 : stateMat(j, k - 1)) - param2 * stateMat(j, k));
						}
					}

					// 3. Go Over two-component interactions
					// A. Interactions in which one of these types takes part, as an input
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l != inter.input1 && l != inter.input2) {
							if (interactionMat(inter.input1,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input1,l)];
								interaction_data &iid = interactionData[interactionMat(inter.input1,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input1) ? iid.conditionalMean2 : iid.conditionalMean1;
								// Handle the case in which the one of the outputs is the other species
								// Check if the other is one of the outputs
								bool isOutput = false;
								for (size_t l = 0; l < ii.outputs.size(); l++) {
									if (ii.outputs[l] == inter.input2) isOutput = true;
								}
								if (isOutput) {
									interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
										((((j == input1.cutoff) || k == 0)? 0 : 
											(j + 1) * means[j+1] * stateMat(j+1, k - 1)) -
										(j * means[j] * stateMat(j,k)));
								} else {
									interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
										(((j == input1.cutoff) ? 0 : 
											(j + 1) * means[j+1] * stateMat(j+1, k)) -
										(j * means[j] * stateMat(j,k)));

								}
							}
							if (interactionMat(inter.input2,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input2,l)];
								interaction_data &iid = interactionData[interactionMat(inter.input2,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input2) ? iid.conditionalMean2 : iid.conditionalMean1;
								bool isOutput = false;
								for (size_t l = 0; l < ii.outputs.size(); l++) {
									if (ii.outputs[l] == inter.input1) isOutput = true;
								}
								if (isOutput) {
									interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
										((((k == input2.cutoff) || j == 0) ? 0 : 
											(k + 1) * means[k+1] * stateMat(j - 1,k + 1)) -
										(k * means[k] * stateMat(j,k)));
								} else {
									interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
										(((k == input2.cutoff) ? 0 : 
											(k + 1) * means[k+1] * stateMat(j, k + 1)) -
										(k * means[k] * stateMat(j,k)));

								}
							}
						}
					}

					// THERE IS A BUG HERE - it may be that an interaction produces both...
					// B. Interactions in which one of these types takes part, as an output
					for (size_t l = 0; l < indexedOutputs[inter.input1].size(); l++) {
						size_t pos = indexedOutputs[inter.input1][l];
						const interaction& ii = interactions[pos];
						const interaction_data &iid = interactionData[pos];
						const species &s1 = chemicalTypes[ii.input1];
						const species &s2 = chemicalTypes[ii.input2];
						// if the other input of this interaction is also
						// an input of the other interaction - we have taken care of this already
						if (ii.input1 != inter.input2 && ii.input2 != inter.input2) {
							// both input species should be traced over, as they are not part of this
							// interaction
							interRes(j,k) += (s1.A + s2.A) * iid.corr *
								((j == 0 ? 0 : stateMat(j - 1, k)) - stateMat(j, k));
						}
					}

					for (size_t l = 0; l < indexedOutputs[inter.input2].size(); l++) {
						size_t pos = indexedOutputs[inter.input2][l];
						const interaction& ii = interactions[pos];
						const interaction_data &iid = interactionData[pos];
						const species &s1 = chemicalTypes[ii.input1];
						const species &s2 = chemicalTypes[ii.input2];
						// if the other input of this interaction is also
						// an input of the other interaction - we have taken care of this already
						if (ii.input1 != inter.input1 && ii.input2 != inter.input1) {
							// both species should be traced over, as they are not part of this
							// interaction
							// CHECK IF IT WAS NOT ALREADY COUNTED
							interRes(j,k) += (s1.A + s2.A) * iid.corr *
								((k == 0 ? 0 : stateMat(j, k - 1)) - stateMat(j, k));
						}

					}

					//4. Go over dissociations
					//A. Dissiciations in which one of the types is the input
					if (dissociationMat(inter.input1, inter.input1) != -1) {
						const dissociation& di = dissociations[dissociationMat(inter.input1, inter.input1)];
						// Check if the other one is an output of the dissociation
						if (dissociationMat(inter.input1, inter.input2) != -1) {
							interRes(j,k) += di.D *
								(((((j + 1) <= input1.cutoff) && k > 0) ? 
									((j+1)*stateMat(j+1, k - 1)) : 0) - (j * stateMat(j,k)));
						} else { 
							interRes(j,k) += di.D *
								(((((j + 1) <= input1.cutoff)) ? 
									((j+1)*stateMat(j+1, k)) : 0) - (j * stateMat(j,k)));
						}
					}
					if (dissociationMat(inter.input2, inter.input2) != -1) {
						const dissociation& di = dissociations[dissociationMat(inter.input2, inter.input2)];
						// Check if the other one is an output of the dissociation
						if (dissociationMat(inter.input2, inter.input1) != -1) {
							interRes(j,k) += di.D *
								(((((k + 1) <= input2.cutoff) && j > 0) ? 
									((k+1)*stateMat(j-1, k+1)) : 0) - (k * stateMat(j,k)));
						} else { 
							interRes(j,k) += di.D *
								(((((k + 1) <= input2.cutoff)) ? 
									((k+1)*stateMat(j, k+1)) : 0) - (k * stateMat(j,k)));
						}
					}
					
					// B. Dissociations with either of these is the output
					// There are several cases - both are the output of the dissociation, one is, or none are					
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l == inter.input1 || l == inter.input2) { 
							// One of the interacting types is the input, continue since we've taken care of this
							continue;
						} 
						int first = dissociationMat(l, inter.input1);
						int second = dissociationMat(l, inter.input2);
						int firstIn = interactionMat(l, inter.input1);
						int secondIn = interactionMat(l, inter.input2);

						double param1 = 0;
						double param2 = 0;
						if (first != -1 || second != -1) {													
							// If the reactants produces only one and does not interact with it, 
							// or produces both and either does not interact with none
							// interactes with both (special case)
							if (((first != -1 && second != -1) && 
								((firstIn != -1 && secondIn != -1) || (firstIn == -1 && secondIn == -1))) ||
								(first != -1 && firstIn == -1) || (second != -1 && secondIn == -1)) {	
									interaction &ii = interactions[dissociations[first].locations[0]];
									interaction_data &iid = interactionData[dissociations[first].locations[0]];
									// Check if it is the first or second input
									param1 = param2 = (l == ii.input1) ?
										iid.mean1 : 
										iid.mean2; 							
							} else if (first != -1 && firstIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input1)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input1)];
								// Check if it is the first or second input								
								if (j != 0) {
									param1 = (l == ii.input1) ?
										iid.conditionalMean1[j - 1] : 
										iid.conditionalMean2[j - 1]; 								
								}
								param2 = (l == ii.input1) ?
									iid.conditionalMean1[j] : 
									iid.conditionalMean2[j]; 							
							} else if (second != -1 && secondIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input2)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input2)];
								// Check if it is the first or second input								
								if (k != 0) {
									param1 = (l == ii.input1) ?
										iid.conditionalMean1[k - 1] : 
										iid.conditionalMean2[k - 1]; 								
								}
								param2 = (l == ii.input1) ?
									iid.conditionalMean1[k] : 
									iid.conditionalMean2[k]; 							

							} else {
								vmerror("BUGBUGBUG");
							}

							const dissociation& di = dissociations[first != -1 ? first : second];
							if (first != -1 && second != -1) { 
								interRes(j,k) += di.D *
									(param1 * (j == 0 || k == 0 ? 0 : stateMat(j - 1, k - 1)) - param2 * stateMat(j, k));
							} else if (first != -1) {
								interRes(j,k) += di.D *
									(param1 * (j == 0 ? 0 : stateMat(j - 1, k)) - param2 * stateMat(j, k));
							} else if (second != -1) { 
								interRes(j,k) += di.D *
									(param1 * (k == 0 ? 0 : stateMat(j, k - 1)) - param2 * stateMat(j, k));
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
	virtual vec compute_derivative_derivative(const vec& state) {
		updateParameters(state);
		vec res = state.copy();

		// First - go over the two-component interactions
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction& inter = interactions[i];
			species& input1 = chemicalTypes[inter.input1];
			species& input2 = chemicalTypes[inter.input2];
			Matrix<double>& interRes = res[i];			

			for (size_t j = 0; j <= input1.cutoff; ++j) {
				for (size_t k = 0; k <= input2.cutoff; ++k) {

					// 1. Take care of the non interaction part and this pair interaction
					interRes(j, k) =
						// Non-interaction parts
						- input1.Flux 
						- input2.Flux
						- input1.W * j 
						- input2.W * k 
						// current pair interaction
						- (input1.A + input2.A) * j * k;

					// 2. Go over self interactions:
					// A. Self interactions where these are the inputs
					if (interactionMat(inter.input1, inter.input1) != -1) {
						interRes(j,k) -= input1.A * j * (j-1);
					}
					
					if (interactionMat(inter.input2, inter.input2) != -1) {
						interRes(j,k) -= input2.A * k * (k-1);
					}

					// B. Self interactions with either of these is the output
					// There are several cases - both are the output of the self interaction, one is, or none are
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l == inter.input1 || l == inter.input2) { 
							// One of the interacting types is the input, continue since we've taken care of this
							continue;
						} 
						int first = selfInteractionMat(l, inter.input1);
						int second = selfInteractionMat(l, inter.input2);
						int firstIn = interactionMat(l, inter.input1);
						int secondIn = interactionMat(l, inter.input2);
						
						double param2 = 0;
						if (first != -1 || second != -1) {													
							// If the reactants produces only one and does not interact with it, 
							// or produces both and either does not interact with none
							// interactes with both (special case)
							if (((first != -1 && second != -1) && 
								((firstIn != -1 && secondIn != -1) || (firstIn == -1 && secondIn == -1))) ||
								(first != -1 && firstIn == -1) || (second != -1 && secondIn == -1)) {	
									interaction &ii = interactions[selfInteractions[first].locations[0]];
									interaction_data &iid = interactionData[selfInteractions[first].locations[0]];
									// Check if it is the first or second input
									param2 = (l == ii.input1) ?
										iid.secondMoment1 - iid.mean1 : 
										iid.secondMoment2 - iid.mean2; 							
							} else if (first != -1 && firstIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input1)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input1)];
								// Check if it is the first or second input																
								param2 = (l == ii.input1) ?
									iid.conditionalSecond1[j] - iid.conditionalMean1[j] : 
									iid.conditionalSecond2[j] - iid.conditionalMean2[j]; 							
							} else if (second != -1 && secondIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input2)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input2)];
								// Check if it is the first or second input								
								param2 = (l == ii.input1) ?
									iid.conditionalSecond1[k] - iid.conditionalMean1[k] : 
									iid.conditionalSecond2[k] - iid.conditionalMean2[k]; 							

							} else {
								vmerror("BUGBUGBUG");
							}
						}
						interRes(j,k) -= chemicalTypes[l].A * param2;
					}

					// 3. Go Over two-component interactions
					// A. Interactions in which one of these types takes part, as an input
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l != inter.input1 && l != inter.input2) {
							if (interactionMat(inter.input1,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input1,l)];
								interaction_data &iid = interactionData[interactionMat(inter.input1,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input1) ? iid.conditionalMean2 : iid.conditionalMean1;
								interRes(j, k) -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
									j * means[j];
							}
							if (interactionMat(inter.input2,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input2,l)];
								interaction_data &iid = interactionData[interactionMat(inter.input2,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input2) ? iid.conditionalMean2 : iid.conditionalMean1;
								interRes(j, k) -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
										k * means[k];
							}
						}
					}

					// THERE IS A BUG HERE - it may be that an interaction produces both...
					// B. Interactions in which one of these types takes part, as an output
					for (size_t l = 0; l < indexedOutputs[inter.input1].size(); l++) {
						size_t pos = indexedOutputs[inter.input1][l];
						const interaction& ii = interactions[pos];
						const interaction_data &iid = interactionData[pos];
						const species &s1 = chemicalTypes[ii.input1];
						const species &s2 = chemicalTypes[ii.input2];
						// if the other input of this interaction is also
						// an input of the other interaction - we have taken care of this already
						if (ii.input1 != inter.input2 && ii.input2 != inter.input2) {
							// both input species should be traced over, as they are not part of this
							// interaction
							interRes(j,k) -= (s1.A + s2.A) * iid.corr;								
						}
					}

					for (size_t l = 0; l < indexedOutputs[inter.input2].size(); l++) {
						size_t pos = indexedOutputs[inter.input2][l];
						const interaction& ii = interactions[pos];
						const interaction_data &iid = interactionData[pos];
						const species &s1 = chemicalTypes[ii.input1];
						const species &s2 = chemicalTypes[ii.input2];
						// if the other input of this interaction is also
						// an input of the other interaction - we have taken care of this already
						if (ii.input1 != inter.input1 && ii.input2 != inter.input1) {
							// both species should be traced over, as they are not part of this
							// interaction
							// CHECK IF IT WAS NOT ALREADY COUNTED
							interRes(j,k) -= (s1.A + s2.A) * iid.corr;								
						}

					}

					//4. Go over dissociations
					//A. Dissiciations in which one of the types is the input
					if (dissociationMat(inter.input1, inter.input1) != -1) {
						const dissociation& di = dissociations[dissociationMat(inter.input1, inter.input1)];
						interRes(j,k) -= di.D * j;
					}
					if (dissociationMat(inter.input2, inter.input2) != -1) {
						const dissociation& di = dissociations[dissociationMat(inter.input2, inter.input2)];
						interRes(j,k) -= di.D * k;
					}

					// B. Dissociations with either of these is the output
					// There are several cases - both are the output of the dissociation, one is, or none are					
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l == inter.input1 || l == inter.input2) { 
							// One of the interacting types is the input, continue since we've taken care of this
							continue;
						} 
						int first = dissociationMat(l, inter.input1);
						int second = dissociationMat(l, inter.input2);
						int firstIn = interactionMat(l, inter.input1);
						int secondIn = interactionMat(l, inter.input2);

						double param2 = 0;
						if (first != -1 || second != -1) {													
							// If the reactants produces only one and does not interact with it, 
							// or produces both and either does not interact with none
							// interactes with both (special case)
							if (((first != -1 && second != -1) && 
								((firstIn != -1 && secondIn != -1) || (firstIn == -1 && secondIn == -1))) ||
								(first != -1 && firstIn == -1) || (second != -1 && secondIn == -1)) {	
									interaction &ii = interactions[dissociations[first].locations[0]];
									interaction_data &iid = interactionData[dissociations[first].locations[0]];
									// Check if it is the first or second input
									param2 = (l == ii.input1) ?
										iid.mean1 : 
										iid.mean2; 							
							} else if (first != -1 && firstIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input1)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input1)];
								param2 = (l == ii.input1) ?
									iid.conditionalMean1[j] : 
									iid.conditionalMean2[j]; 							
							} else if (second != -1 && secondIn != -1) {
								interaction &ii = interactions[interactionMat(l, inter.input2)];
								interaction_data &iid = interactionData[interactionMat(l, inter.input2)];
								// Check if it is the first or second input								
								param2 = (l == ii.input1) ?
									iid.conditionalMean1[k] : 
									iid.conditionalMean2[k]; 							
							} else {
								vmerror("BUGBUGBUG");
							}

							const dissociation& di = dissociations[first != -1 ? first : second];
							interRes(j,k) -= di.D * param2;
						}
					}
				}
			}
		}
		return res;
	}

	/** 
    	 * Normalize the results (if relevant)
	 * 
	 * @param state The un-normalized state
	 * @return normalized result
	 */ 
	virtual vec normalize(const vec& state) { 
		vec res = state.copy();
		for (size_t i = 0; i < interactions.size(); i++) { 
			double sum = 0.0;
			for (size_t n = 0; n < state[i].nrows(); n++) {
				sum += state[i].row(n).sum();
			}
			res[i] /= sum;
		}
		return res;
	}

	vec createInitialConditions() {
		vec result(interactions.size());
		// Start with empty state;
		for (size_t i = 0; i < interactions.size(); i++) {
			result[i].resize(
				chemicalTypes[interactions[i].input1].cutoff + 1,
				chemicalTypes[interactions[i].input2].cutoff + 1);

			result[i] = 0.0;
			result[i](0,0) = 1.0;
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
			const interaction_data &id = interactionData[i];
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
