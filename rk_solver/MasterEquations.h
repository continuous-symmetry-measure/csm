
#ifndef MASTER_EQ_H
#define MASTER_EQ_H

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
	size_t cutoff;			// Cutoff (max number of particles)
};

struct interaction_desc {
	string input2;			// The first species' name
	string input1;			// The second species' name
	string output;			// The output material	
};

class MasterEquations : public EquationSet<double>, public RKResultProcessor<double>{
private:
	
	struct interaction {
		size_t input1;			// The first input type
		size_t input2;			// The second input type
		size_t output;			// The output type 		
	};

	struct self_interaction {
		size_t output;			// The output of this interaction
		size_t input;			// The input of this interaction
		vector<size_t> locations;	// The interactions containing this self interaction
	};

	typedef map<string, int> name_index_map;

	name_index_map indexer;				// A converter from a species' name to its index in the 
							// chemicalTypes and meanValueHelper vectors
	name_index_map no_inter_map;			// A map of the non-interacting types
	
	vector<interaction> interactions;	// The interactions 
	vector<species> chemicalTypes;				// The chemical types	
	vector<self_interaction> selfInteractions;	// The self interactions
	Matrix<int> interactionMat;			// A matrix of interactions
	vector<string> outputTypes;

	// Assume output can only come from one source !!!
	Vector<int> indexedOutputs;
	Vector<int> indexedSelfOutputs;	

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

	string getOutputName(int index) { 
		return ((size_t)index >= chemicalTypes.size()) ? outputTypes[index - chemicalTypes.size()] : chemicalTypes[index].name;
	}

public:
	virtual ~MasterEquations() {}
	MasterEquations(const vector<species>& types, const vector<interaction_desc>& unprocessedInteractions) {

		indexedOutputs.resize(types.size());
		indexedOutputs = -1;
		indexedSelfOutputs.resize(types.size());
		indexedSelfOutputs = -1;
		interactionMat.resize(types.size(), types.size());
		interactionMat = -1;	

		// First check up on all the types of species		
		chemicalTypes = types;
		int pos = 0;
		maxPos = 1;
		
		for (size_t i = 0; i < types.size(); i++) {
			if (indexer.find(types[i].name) != indexer.end()) {
				cerr << "Species " << types[i].name << " appears more than once";
				exit(1);
			} else {
				indexer[types[i].name] = pos;
				pos++;
			}
			maxPos *= (types[i].cutoff + 1);
			if (i == 0) {		
				diffs.push_back(1);
			} else {
				diffs.push_back(diffs[i - 1] * (types[i - 1].cutoff + 1));
			}
			cout << diffs[i] << " ";
		}
		cout << endl;

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

			bool isInteracting = false;
			int outputIndex = 0;
			// The output type can be a non-interacting type
			if (outputType == indexer.end()) {
				name_index_map::iterator output = no_inter_map.find(input.output);
				if (output == no_inter_map.end()) { 
					outputIndex = chemicalTypes.size() + no_inter_map.size();
					no_inter_map[input.output] = outputIndex;
					outputTypes.push_back(input.output);
				} else {	
					outputIndex = output->second;
				}
			} else {
				isInteracting = true;
				outputIndex = outputType->second;
			}

			// Everything is OK - fill parameters
			if (first->second == second->second) {
				// This is a self interaction
				self_interaction si;
				si.input = first->second;
				si.output = outputIndex;
				selfInteractions.push_back(si);				
				interactionMat(first->second, first->second) = selfInteractions.size() - 1;
				if (isInteracting) {
					indexedSelfOutputs[si.output] = selfInteractions.size() - 1;
				}
			} else {
				interaction output;		
				output.input1 = first->second;
				output.input2 = second->second;
				output.output = outputIndex;
				interactions.push_back(output);
				interactionMat(output.input1, output.input2) = interactions.size() - 1;
				interactionMat(output.input2, output.input1) = interactions.size() - 1;
				if (isInteracting) {
					indexedOutputs[output.output] = interactions.size() - 1;
				}
			}
		}

		// Now, go over the self interactions, and index the interactions containing them
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];			
			for (size_t j = 0; j < chemicalTypes.size(); j++) { 
				if (si.input != j && interactionMat(si.input, j) != -1) { 
					si.locations.push_back(interactionMat(si.input, j));
				}	
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
	virtual void updateParameters(const vec& state) {}

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
					// Diffusion Term
					chemicalTypes[i].Diffusion * 
					(((indices[i] == chemicalTypes[i].cutoff) ? 0 : (indices[i] + 1) * state[pos + diffs[i]]) - 
					(indices[i] * state[pos]));
			}			


			for (size_t i = 0; i < selfInteractions.size(); i++) {
				const self_interaction &si = selfInteractions[i];
				size_t ii = indices[si.input];
				size_t oi = indices[si.output];
				if (si.output >= chemicalTypes.size()) { 
					// The output of this interaction is a non-interacting species
					res[pos] += chemicalTypes[si.input].A *
					    (((ii >= chemicalTypes[si.input].cutoff - 1) ? 0 : 
						((ii + 2) * (ii + 1) * state[pos + (2 * diffs[si.input])])) -
					    (ii * (ii - 1) * state[pos]));				
				} else {	
					res[pos] += chemicalTypes[si.input].A *
					    ((((ii >= chemicalTypes[si.input].cutoff - 1) || 
					      (oi == 0)) ? 0 : ((ii + 2) * (ii + 1) * state[pos + (2 * diffs[si.input]) - diffs[si.output]])) - 
					    (ii * (ii - 1) * state[pos]));				
				}				
			}

			// Next - Go over two-component interactions
			for (size_t i = 0; i < interactions.size(); i++) { 
				const interaction &ii = interactions[i];
				size_t i1 = indices[ii.input1];
				size_t i2 = indices[ii.input2];
				size_t oi = indices[ii.output];
				if (ii.output >= chemicalTypes.size())  { 
					// The output of this interaction is a non-interacting species
					res[pos] += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
						((((i1 == chemicalTypes[ii.input1].cutoff) || (i2 == chemicalTypes[ii.input2].cutoff)) ? 
						 	0 : ((i1 + 1) * (i2 + 1) * state[pos + diffs[ii.input1] + diffs[ii.input2] ])) - 
						(i1 * i2 * state[pos]));
				} else { 
					res[pos] += (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
						((((i1 == chemicalTypes[ii.input1].cutoff) || (i2 == chemicalTypes[ii.input2].cutoff) || (oi == 0)) ? 
						 	0 : ((i1 + 1) * (i2 + 1) * state[pos + diffs[ii.input1] + diffs[ii.input2] - diffs[ii.output]])) - 
						(i1 * i2 * state[pos]));
				}
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
					// Diffusion Term
					chemicalTypes[i].Diffusion * indices[i];
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
	static MasterEquations* parseMasterEquations(istream& input) {
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

		return new MasterEquations(types, inters);

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
			file << "Production Rate of " << getOutputName(si.output) << " is " 
				<< (chemicalTypes[si.input].A * (momentVec[si.input] - avgVec[si.input])) << endl;
		}

		for (size_t i = 0; i < interactions.size(); i++) { 
			interaction &ii = interactions[i];	
			file << "Production Rate of " << getOutputName(ii.output) << " is " 
				<< ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input1].A) * corrVec[i])
				<< endl;
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

		file << "Solving Started" << endl;
		for (size_t i = 0; i < chemicalTypes.size(); i++) { 
			file << "<" << chemicalTypes[i].name << ">" << "\t\t";
		}

		for (size_t i = 0; i < interactions.size(); i++) { 
			interaction &ii = interactions[i];	
			file << "R_" << getOutputName(ii.output) << "\t";
		}

		for (size_t i = 0; i < selfInteractions.size(); i++) { 
			self_interaction &si = selfInteractions[i];
			file << "R_" << getOutputName(si.output) << "\t";
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
			file << (chemicalTypes[si.input].A * (momentVec[si.input] - avgVec[si.input]))
				<< "\t";
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
