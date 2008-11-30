
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

class ChemicalNetwork : public EquationSet<double>, public RKResultProcessor<double>{
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

	string getOutputName(int index) { 
		return ((size_t)index >= chemicalTypes.size()) ? outputTypes[index - chemicalTypes.size()] : chemicalTypes[index].name;
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
		for (int i = 0; i < chemicalTypes.size(); i++) { 
			const species& s = chemicalTypes[i];
			res[i] += (s.F - s.W * state[i]);
		
			// Check if there is a self interaction for this:
			if (interactionMat(i,i) != -1) { 
				res[i] -= 2 * s.A * 
				(state[chemicalTypes.size() + interactions.size() + interactionMat(i,i)] - state[i]);
			}
			
			// Check all interactions in which this is the input
			for (int j = 0; j < chemicalTypes.size()) { 
				if (i != j && interactionMat(i,j) != -1) { 
					int pos = interactionMat(i,j);
					const interaction &ii = interactions[pos];
					res[i] -= (chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
							state(chemicalTypes.size() + pos);
				}
			}

			// Check all self interactions in which this is output (currently - only one...)
			if (indexedSelfOutputs[i] != -1) { 
				int pos = indexedSelfOutputs[i];
				const self_interaction &si = selfInteractions[pos];
				res[i] += chemicalTypes[si.input].A * 
					(state[chemicalTypes.size() + interactions.size() + pos] - state[i]);
			}	

			// Check all  interactions in which this is output (currently - only one...)
			if (indexedOutputs[i] != -1) { 
				int pos = indexedOutputs[i];
				interaction &ii = interactions[pos];
				res[i] += ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) *
					(state[chemicalTypes.size() + pos]));
			}				
		}

		// Go over all interactions
		for (int i = 0; i < interactions.size(); i++) { 
			const interaction& ii = interactions[i];
			const species &s1 = chemicalTypes[ii.input1];
			const species &s2 = chemicalTypes[ii.input2];
			res[chemicalTypes.size() + i] = 
				(s1.F * state[ii.input2] + s2.F * state[ii.input1] - 
				 (s1.W + s2.W + s1.A + s1.A) * state[chemicalTypes.size() + i]);
		}
		
		for (int i = 0; i < selfInteractions.size(); i++) { 
			const self_interaction &si = selfInteractions[i];	
			const species &s = chemicalTypes[ii.input];
			res[chemicalTypes.size() + interactions.size() + i] = 
				s.F + 2 * s.F * state[si.input] + s.W * state[si.input] - 
				2 * s.W	* state[chemicalTypes.size() + interactions.size() + i] -
				4 * s.A * (state[chemicalTypes.size() + interactions.size() + i] - state[si.input]);
				
			
			// Check all interactions in which this is the input
			for (int j = 0; j < chemicalTypes.size()) { 
				if (si.input != j && interactionMat(si.input,j) != -1) { 
					int pos = interactionMat(si.input, j);
					const interaction &ii = interactions[pos];
					res[chemicalTypes.size() + interactions.size() + i] -= 
						(chemicalTypes[ii.input1].A + chemicalTypes[ii.input2].A) * 
							state(chemicalTypes.size() + pos);
				}
			}

			// Check all self interactions in which this is output (currently - only one...)
			if (indexedSelfOutputs[si.input] != -1) { 
				int pos = indexedSelfOutputs[si.input];
				const self_interaction &si2 = selfInteractions[pos];
				res[chemicalTypes.size() + interactions.size() + i] += 
					chemicalTypes[si2.input].A * 
					(state[chemicalTypes.size() + interactions.size() + pos] - state[si2.input]);
			}	

			// Check all  interactions in which this is output (currently - only one...)
			if (indexedOutputs[i] != -1) { 
				int pos = indexedOutputs[i];
				interaction &ii = interactions[pos];
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
		file.open("master.out");
		file << "Solving Started" << endl;
	}

	/** 
	 * A single complete RK step has been performed
 	 * @param time	The time
	 * @param dt	The chosen delta
	 * @param state The state after the step
	 */
	virtual void stepPerformed(double time, double dt, const vec& state, const vec& prevState) { 			

		file << "Time: " << time << endl;	
		for (size_t i = 0; i < chemicalTypes.size(); i++) { 
			file << "Mean of " << chemicalTypes[i].name << " is " << state[i] << endl;
		}
		for (size_t i = 0; i < interactions.size(); i++) { 
			interaction &ii = interactions[i];	
			file << "Production Rate of " << getOutputName(ii.output) << " is " 
				<< ((chemicalTypes[ii.input1].A + chemicalTypes[ii.input1].A) * 
					state[i + chemicalTypes.size()])
				<< endl;
		}
		for (size_t i = 0; i < selfInteractions.size(); i++) { 
			self_interaction &si = selfInteractions[i];
			file << "Production Rate of " << getOutputName(si.output) << " is " 
				<< (chemicalTypes[si.input].A * 
					(state[i + chemicalTypes.size() + interactions.size()] - state[si.input]) 
				<< endl;
		}


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
