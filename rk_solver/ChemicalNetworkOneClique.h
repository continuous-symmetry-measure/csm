
#ifndef CHEM_NET
#define CHEM_NET

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>

// Network file format
#define BEGIN_FILE_LABEL	"Begin"
#define END_FILE_LABEL		"End"
#define SPECIES_LABEL		"Species"
#define FLUX_LABEL			"Flux"
#define DIFFUSION_LABEL		"Diffusion"
#define SWEEP_RATE_LABEL	"SweepRate"
#define CUTOFF_LABEL		"Cutoff"
#define INTERACTION_LABEL	"Interaction"
#define ARROW_LABEL			"=>"

// TODO what about self interactions's output
// TODO feedback loops 
// TODO more than one output types

using namespace std;

struct species {
	string name;			// The species name
	double Flux;			// Incoming Flux
	double Diffusion;		// Diffusion Rate	
	double A;			// Sweeping Rate
	int    cutoff;			// Cutoff (max number of particles)
};

struct interaction_desc {
	string input2;			// The first species' name
	string input1;			// The second species' name
	string output;			// The output material	
};

// Here the state is a vector of matrices, each (cutoff1 + 1) x (cutoff2 + 1) 
// sized according to the interaction types
class ChemicalNetwork : public EquationSet<Matrix<double> >, public RKResultProcessor<Matrix<double> >{
private:
	struct interaction {
		size_t input1;			// The first input type
		size_t input2;			// The second input type
		size_t output;			// The output type 		
		Vector<double> mean1;		// The mean values for the first input given the second value
		Vector<double> mean2;		// The mean values for the second input given the first value	
		double totalMean1;		// The total mean of the first input
		double totalMean2;		// The total mean of the second input	
		double corr;			// <input1 * input2>
	};

	struct self_interaction {
		size_t output;			// The output of this interaction
		size_t input;			// The input of this interaction
		double mean;			// The First moment of the input
		double secondMoment;		// The second moment of the input
	};

	typedef map<string, int> name_index_map;

	name_index_map indexer;				// A converter from a species' name to its index in the 
							// chemicalTypes and meanValueHelper vectors
	
	vector<interaction> interactions;	// The interactions 
	vector<species> chemicalTypes;				// The chemical types	
	vector<self_interaction> selfInteractions;	// The self interactions
	Vector<Vector<double> > meanValueHelpers;	// vectors of [0..cutoff] for each species		
	Matrix<int> interactionMat;			// A matrix of interactions

	// Assume output can only come from one source !!!
	Vector<int> indexedOutputs;
	Vector<int> indexedSelfOutputs;	

	ofstream file;

	void print(const interaction& inter) const {	
		std::cout << chemicalTypes[inter.input1].name << " " << chemicalTypes[inter.input2].name << " => " << chemicalTypes[inter.output].name << std::endl;
	}

	void print(const self_interaction& inter) const {	
		std::cout << chemicalTypes[inter.input].name << " " << chemicalTypes[inter.input].name << " => " << chemicalTypes[inter.output].name << std::endl;
	}
	
	Vector<double> prepareResultsVec() { 
		Vector<double> results(2 * selfInteractions.size() + 3 * interactions.size());
		int pos = 0;
		for (size_t i = 0; i < interactions.size(); i++) {
			interaction &inter = interactions[i];
			results[pos++] = inter.totalMean1;
			results[pos++] = inter.totalMean2;
			results[pos++] = inter.corr;
		}
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &inter = selfInteractions[i];
			results[pos++] = inter.mean;
			results[pos++] = inter.secondMoment;
	
		}
		return results;
	}

public:
	virtual ~ChemicalNetwork() {}
	ChemicalNetwork(const vector<species>& types, const vector<interaction_desc>& unprocessedInteractions) {

		indexedOutputs.resize(types.size());
		indexedOutputs = -1;
		indexedSelfOutputs.resize(types.size());
		indexedSelfOutputs = -1;
		interactionMat.resize(types.size(), types.size());
		interactionMat = -1;	

		// First check up on all the types of species		
		chemicalTypes = types;
		meanValueHelpers.resize(types.size());
		int pos = 0;
		for (size_t i = 0; i < types.size(); i++) {
			if (indexer.find(types[i].name) != indexer.end()) {
				cerr << "Species " << types[i].name << " appears more than once";
				exit(1);
			} else {
				indexer[types[i].name] = pos;
				pos++;
			}
			meanValueHelpers[i].resize(types[i].cutoff + 1);
			for (size_t j = 0; j < meanValueHelpers[i].size(); ++j) {
				meanValueHelpers[i][j] = (double)j;
			}
		}
		// Index the interactions as well
		for (size_t i = 0; i < unprocessedInteractions.size(); i++) {
			const interaction_desc& input = unprocessedInteractions[i];
				
			name_index_map::iterator first = indexer.find(input.input1);
			name_index_map::iterator second = indexer.find(input.input2);
			name_index_map::iterator outputType = indexer.find(input.output);
			if (first == indexer.end()) {
				cerr << "Species " << input.input1 << " is not defined";
				exit(1);
			}
			if (second == indexer.end()) {
				cerr << "Species " << input.input2 << " is not defined";
				exit(1);
			}
			if (outputType == indexer.end()) {
				cerr << "Species " << input.output << " is not defined";
				exit(1);
			}

			// Everything is OK - fill parameters
			if (first->second == second->second) {
				// This is a self interaction
				self_interaction si;
				si.input = first->second;				si.output = outputType->second;
				selfInteractions.push_back(si);				
				interactionMat(first->second, first->second) = selfInteractions.size() - 1;
				indexedSelfOutputs[si.output] = selfInteractions.size() - 1;
			} else {
				interaction output;		
				output.input1 = first->second;
				output.input2 = second->second;
				output.output = outputType->second;				
				output.mean1.resize(chemicalTypes[output.input2].cutoff + 1);
				output.mean2.resize(chemicalTypes[output.input1].cutoff + 1);
				interactions.push_back(output);
				interactionMat(output.input1, output.input2) = interactions.size() - 1;
				interactionMat(output.input2, output.input1) = interactions.size() - 1;
				indexedOutputs[output.output] = interactions.size() - 1;
			}
		}
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
			const Matrix<double> &probs = state[i];		

			// Compute means
			inter.totalMean1 = 0.0;
			inter.totalMean2 = 0.0;
			inter.mean1   = 0.0;
			inter.mean2   = 0.0;
			inter.corr    = 0.0;
			double rowSum = 0.0;
			double colSum = 0.0;

			for (int j = 0; j <= chemicalTypes[inter.input1].cutoff; j++) { 
				rowSum = probs.row(j).sum();
				for (int k = 0; k <= chemicalTypes[inter.input2].cutoff; k++) {
					colSum = probs.column(k).sum();
					inter.mean1[k] += (colSum > 0.0) ? (j * probs(j,k) / colSum) : 0.0;	
					inter.mean2[j] += (rowSum > 0.0) ? (k * probs(j,k) / rowSum) : 0.0;
					inter.totalMean1 += j * probs(j,k);
					inter.totalMean2 += k * probs(j,k);					
					inter.corr += k * j * probs(j,k);
				}				
			}	

			// Sum over all just to test...
			double sum = 0;
			for (size_t j = 0; j < probs.nrows(); j++) {
				sum += probs.row(j).sum();
			}
			cout << endl << "Interaction #" << (i + 1) << " sum: " << sum <<  endl;
		}
	
		// Go over the self interactions
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction& inter = selfInteractions[i];
			const Matrix<double> &probs = state[interactions.size() + i];
			inter.mean = DotProduct(probs.row(0), meanValueHelpers[inter.input]);		
		inter.secondMoment = DotProduct(probs.row(0), (meanValueHelpers[inter.input] * meanValueHelpers[inter.input]));
			cout << endl << "Self Interaction #" << (i+1) << " sum: " << probs.row(0).sum() << endl;
		}
	}

	/** 
 	 * Compute the error, given the previous state, the current state, and the used delta
  	 */
	virtual double computeError(const vec& state, const vec& prevState, double usedDelta) { 
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
		vec res = state.copy();		
			
		// First - go over the two-component interactions
		for (size_t i = 0; i < interactions.size(); i++) {	
			interaction& inter = interactions[i];
			species& input1 = chemicalTypes[inter.input1];
			species& input2 = chemicalTypes[inter.input2];
			Matrix<double>& interRes = res[i];
			
			const Matrix<double>& stateMat = state[i];
			for (int j = 0; j <= input1.cutoff; ++j) {
				for (int k = 0; k <= input2.cutoff; ++k) {

					// 1. Take care of the non interaction part and this pair interaction					
					interRes(j, k) = 
						// Non-interaction parts
						input1.Flux * ((j == 0 ? 0 : stateMat(j - 1, k))  - stateMat(j, k)) + 
						input2.Flux * ((k == 0 ? 0 : stateMat(j, k - 1))  - stateMat(j, k)) + 
						input1.Diffusion * (((j == input1.cutoff)? 0 : ((j + 1) * stateMat(j + 1, k))) - (j * stateMat(j, k))) + 
						input2.Diffusion * (((k == input2.cutoff)? 0 : ((k + 1) * stateMat(j, k + 1))) - (k * stateMat(j, k))) + 
						// current pair interaction
						(input1.A + input2.A) * (((j == input1.cutoff) || (k == input2.cutoff) ? 0 : 
							((j+1)*(k+1)*stateMat(j+1, k+1))) - (j * k * stateMat(j,k)));					

					// 2. Go over self interactions:
					// A. Self interactions where these are the inputs
					if (interactionMat(inter.input1, inter.input1) != -1) {			
						interRes(j,k) += input1.A * 
							((((j + 2) <= input1.cutoff) ? ((j+1)*(j+2)*stateMat(j+2, k)) : 0) - (j * (j-1)*stateMat(j,k)));
					}
					if (interactionMat(inter.input2, inter.input2) != -1) {					
						interRes(j,k) += input2.A * 
							((((k + 2) <= input2.cutoff) ? ((k+1)*(k+2)*stateMat(j, k+2)) : 0) - (k * (k-1)*stateMat(j,k)));
					}

					// B. Self interactions with either of these is the output
					if (indexedSelfOutputs[inter.input1] != -1) { 
						int pos = indexedSelfOutputs[inter.input1];
						self_interaction& si = selfInteractions[pos];
						interRes(j,k) += (si.secondMoment - si.mean) * chemicalTypes[si.input].A * 
							((j == 0 ? 0 : stateMat(j - 1, k)) - stateMat(j, k));
								
					}
					if (indexedSelfOutputs[inter.input2] != -1) { 
						int pos = indexedSelfOutputs[inter.input2];
						self_interaction& si = selfInteractions[pos];
						interRes(j,k) += (si.secondMoment - si.mean) * chemicalTypes[si.input].A * 
							((k == 0 ? 0 : stateMat(j, k - 1)) - stateMat(j, k));						
					} 		
	 

					// 3. Go Over two-component interactions				
					// A. Interactions in which one of these types takes part, as an input
					for (size_t l = 0; l < chemicalTypes.size(); l++) {
						if (l != inter.input1 && l != inter.input2) {							
							if (interactionMat(inter.input1,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input1,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input1) ? ii.mean2 : ii.mean1;
								interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
									(((j == input1.cutoff) ? 0 : (j + 1) * means[j+1] * stateMat(j+1, k)) - 
										(j * means[j] * stateMat(j,k)));
							} 
							if (interactionMat(inter.input2,l) != -1) {
								interaction &ii = interactions[interactionMat(inter.input2,l)];
								// Check if it is the first or second input
								Vector<double> &means = (ii.input1 == inter.input2) ? ii.mean2 : ii.mean1;
								interRes(j, k) += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
									(((k == input2.cutoff) ? 0 : (k + 1) * means[k+1] * stateMat(j, k + 1)) - 
									(k * means[k] * stateMat(j,k)));

							}  
						}
					}

					// B. Interactions in which one of these types takes part, as an output (there can be only one)
					if (indexedOutputs[inter.input1] != -1) { 
						interaction& ii = interactions[indexedOutputs[inter.input1]];
						species &s1 = chemicalTypes[ii.input1];
						species &s2 = chemicalTypes[ii.input2];
						// There are two cases - either one of the input species is the second one, 
						// or not						
						if (ii.input1 == inter.input2 || ii.input2 == inter.input2) { 	
							// if it is, than only the third species needs to 
							// be traced over
							double tracedMean = (ii.input1 == inter.input2) ? ii.totalMean2 : ii.totalMean1;
							interRes(j,k) += (s1.A + chemicalTypes[ii.input2].A) * 
								tracedMean * 
								(((j == 0 || (k == chemicalTypes[inter.input2].cutoff)) ? 0 : 
									((k + 1) * stateMat(j - 1, k + 1))) - k * stateMat(j, k));
						} else { 
							// otherwise - both species should be traced over, as they are not part of this 
							// interaction
							interRes(j,k) += (s1.A + s2.A) * ii.corr * 
								((k == 0 ? 0 : stateMat(j, k - 1)) - stateMat(j, k));
						}						
								
					}

					if (indexedOutputs[inter.input2] != -1) { 
						interaction& ii = interactions[indexedOutputs[inter.input2]];
						species &s1 = chemicalTypes[ii.input1];
						species &s2 = chemicalTypes[ii.input2];
						// There are two cases - either one of the input species is the second one, 
						// or not						
						if (ii.input1 == inter.input1 || ii.input2 == inter.input1) { 	
							// if it is, than only the third species needs to 
							// be traced over
							double tracedMean = (ii.input1 == inter.input1) ? ii.totalMean2 : ii.totalMean1;
							interRes(j,k) += (s1.A + s2.A) * 
								tracedMean * 
								(((k == 0 || (j == chemicalTypes[inter.input1].cutoff)) ? 0 : 
									((j + 1) *  stateMat(j + 1, k - 1))) - j * stateMat(j, k));
						} else { 
							// otherwise - both species should be traced over, as they are not part of this 
							// interaction
							interRes(j,k) += (s1.A + s2.A) * ii.corr *
								((j == 0 ? 0 : stateMat(j - 1, k)) - stateMat(j, k));
						}						
								
					}						

				}
			}
		}

		// Next - go over the self interactions
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			// Read Parameters
			Matrix<double> interRes = res[interactions.size() + i];	
			Matrix<double> stateMat = state[interactions.size() + i];	
			self_interaction& inter = selfInteractions[i];
			species &input = chemicalTypes[inter.input];
				
			// 1. Take care of the non-interaction parts and this self-interaction 
			for (int j = 0; j <= input.cutoff; j++) {
				
				interRes(0, j) = input.Flux * ((j == 0 ? 0 : stateMat(0, j - 1))  - stateMat(0, j)) +
					input.Diffusion * (((j == input.cutoff)? 0 : ((j + 1) * stateMat(0, j + 1))) - (j * stateMat(0, j))) +		
					input.A * ((((j + 2) <= input.cutoff) ? ((j+1)*(j+2)*stateMat(0, j+2)) : 0) - (j * (j-1)*stateMat(0,j)));
			}
		
			// 2. Go over all interactions where this may be the input
			for (size_t k = 0; k < chemicalTypes.size(); k++) {	
				if (k != inter.input && interactionMat(inter.input, k) != -1) {
					species &s = chemicalTypes[k];
					interaction &ii = interactions[interactionMat(inter.input, k)];	
					double mean = 0.0;
					if (ii.input1 == k) {
						mean = ii.totalMean1;
					} else {
						mean = ii.totalMean2;
					}
					for (int j = 0; j <= input.cutoff; j++) { 
						interRes(0, j) += (input.A + s.A) * mean *
								(((j == input.cutoff) ? 0 : (j + 1) * stateMat(0, j + 1)) - 
									(j * stateMat(0,j)));	
					} 						
				}
			}
				
			// 3. Go over the interaction where this may be the output
			if (indexedOutputs[inter.input] != -1) {
				interaction &ii = interactions[indexedOutputs[inter.input]];	
				species &s1 = chemicalTypes[ii.input1];
				species &s2 = chemicalTypes[ii.input2];
				for (int j = 0; j <= input.cutoff; j++) { 					
					interRes(0, j) += (s1.A + s2.A) * ii.corr * 
						(((j == 0) ? 0 : stateMat(0, j - 1)) - stateMat(0,j));	
				}
			}
		
			// 4. Go over all self interactions where this may be the output
			if (indexedSelfOutputs[inter.input] != -1) {
				self_interaction &ii = selfInteractions[indexedSelfOutputs[inter.input]];	
				species &s = chemicalTypes[ii.input];
				for (int j = 0; j <= input.cutoff; j++) { 				
					interRes(0, j) += (s.A) * (ii.secondMoment - ii.mean) * 
						(((j == 0) ? 0 : stateMat(0, j - 1)) - stateMat(0,j));						
				} 
			}
		}
		return res;
	}

	vec createInitialConditions() {
		vec result(interactions.size() + selfInteractions.size());
		// Start with empty state;
		for (size_t i = 0; i < interactions.size(); i++) {
			result[i].resize(
				chemicalTypes[interactions[i].input1].cutoff + 1,
				chemicalTypes[interactions[i].input2].cutoff + 1);
		
			result[i] = 0.0;
			result[i](0,0) = 1.0;		
		}
		for (size_t i=0; i < selfInteractions.size(); i++) {
			result[interactions.size() + i].resize(1, chemicalTypes[selfInteractions[i].input].cutoff + 1);
			result[interactions.size() + i] = 0.0;
			result[interactions.size() + i](0,0) = 1.0;
		}
		return result;
	}

	/** 
	 * Read and parse a chemical network from an input stream. 
	 * The file structure is as follows:
	 * 
	 * Begin
	 * Species		Flux		Diffusion		SweepRate	Cutoff
	 * <name>		0.4		8			4		6
	 * ...
	 * 
	 * Interaction
	 * <input1> <input2> => <output1> 
	 * <input1> <input2> => <output1> 
	 * End
	 *
	 * @param input The input stream
	 * @return a newly-created chemical network
	 */
	static ChemicalNetwork* parseChemicalNetwork(istream& input) {
		vector<species> types;
		vector<interaction_desc> inters;
		string str1, str2, str3, str4, str5;		
		input >> skipws;
		input >> str1;
		if (str1 != BEGIN_FILE_LABEL) {
			cerr << "File should begin with 'Begin'" << endl;
			return NULL;
		}
		input >> str1 >> str2 >> str3 >> str4 >> str5;
		if (str1 != SPECIES_LABEL || str2 != FLUX_LABEL || str3 != DIFFUSION_LABEL || str4 != SWEEP_RATE_LABEL || str5 != CUTOFF_LABEL) {
			cerr << "File does not contain 'Species Flux Diffusion SweepRate Cutoff' header" << endl;
			return NULL;
		}
		input >> str1;
		while (str1 != INTERACTION_LABEL) {	
			species sp;
			sp.name = str1;
			input >> sp.Flux >> sp.Diffusion >> sp.A >> sp.cutoff;
			types.push_back(sp);
			input >> str1;
			cerr << str1 << endl;
		}
		input >> str1;
		while (str1 != END_FILE_LABEL) {
			interaction_desc inter;
			input >> str2 >> str3 >> str4;
			if (str3 != ARROW_LABEL) {
				cerr << "Interaction format: <input1> <input2> => <output1>" << endl;
				return NULL;
			}
			inter.input1 = str1;
			inter.input2 = str2;
			inter.output = str4;			
			inters.push_back(inter);
			input >> str1;
					cerr << "analyzing interaction" << endl;
		}

		return new ChemicalNetwork(types, inters);

	}

	/**
         * Announce that the solving has started
	 */
	virtual void solvingStarted(const rk_params &params, const vec& initialState) { 
		file.open("pairs.out");
		file << "Solving Started" << endl;
	}

	/** 
	 * A single complete RK step has been performed
 	 * @param time	The time
	 * @param dt	The chosen delta
	 * @param state The state after the step
	 */
	virtual void stepPerformed(double time, double dt, const vec& state, const vec& prevState) { 			
		file << "Time: " << time << endl ;	
		for (size_t i = 0; i < interactions.size(); ++i) {	
			interaction& inter = interactions[i];
			file << "Mean of " << chemicalTypes[inter.input1].name << " is " << inter.totalMean1 << endl;
			file << "Mean of " << chemicalTypes[inter.input2].name << " is " << inter.totalMean2 << endl;
			file << "ProductionRate of " << chemicalTypes[inter.output].name << " is " 
				<< (chemicalTypes[inter.input1].A + chemicalTypes[inter.input2].A)*(inter.corr) << endl;
		}		
		for (size_t i = 0; i < selfInteractions.size(); ++i) {	
			self_interaction& inter = selfInteractions[i];
			file << "Mean of " << chemicalTypes[inter.input].name << " is " << inter.mean << endl;
			file << "ProductionRate of " << chemicalTypes[inter.output].name << " is " 
				<< chemicalTypes[inter.input].A * (inter.secondMoment - inter.mean) << endl;
		}		
//		file << time << " " << interactions[0].totalMean1 << " "  << interactions[0].totalMean2 << endl;

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
