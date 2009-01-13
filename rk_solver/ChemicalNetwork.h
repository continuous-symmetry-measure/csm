#ifndef CHEM_NET
#define CHEM_NET

#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <vm/vec_mat.h>

// Network file format
#define BEGIN_FILE_LABEL	"Begin"
#define END_FILE_LABEL		"End"
#define SPECIES_LABEL		"Species"
#define FLUX_LABEL			"Flux"
#define DESORPTION_LABEL	"Desorption"
#define SWEEP_RATE_LABEL	"SweepRate"
#define CUTOFF_LABEL		"Cutoff"
#define INTERACTION_LABEL	"Interaction"
#define ARROW_LABEL			"=>"
#define DISSOCIATION_LABEL	"Dissociation"
#define SEPARATOR_LABEL		"|"

// TODO what about self interactions's output
// TODO feedback loops
// TODO more than one output types
using namespace std;

struct species {
	string name;			// The species name
	double Flux;			// Incoming Flux
	double W;				// W Rate
	double A;				// Sweeping Rate
	size_t cutoff;			// Cutoff (max number of particles)
};

struct interaction_desc {
	string input2;			// The first species' name
	string input1;			// The second species' name
	vector<string> outputs;			// The output material
};

struct dissociation_desc {
	string input;
	vector<string> outputs;
	double D;
};

struct parsed_network {
	vector<species> types;
	vector<dissociation_desc> dis;
	vector<interaction_desc> interactions;
};


class ChemicalNetwork {
protected:
	struct interaction {
		size_t input1;			// The first input type
		size_t input2;			// The second input type
		vector<size_t> outputs;			// The output type
	};

	struct self_interaction {
		size_t input;					// The input of this interaction
		vector<size_t> outputs;			// The output of this interaction
		vector<size_t> locations;		// The interactions containing this self interaction
	};

	struct dissociation {
		size_t input;					// The input type
		vector<size_t> outputs;			// The outputs
		double D;						// The dissociation rate
		vector<size_t> locations;		// The interactions containing this self interaction
	};

	typedef map<string, int> name_index_map;

	name_index_map indexer;				// A converter from a species' name to its index in the
										// chemicalTypes and meanValueHelper vectors
	name_index_map no_inter_map;		// A map of the non-interacting types


	vector<species> chemicalTypes;				// The chemical types
	vector<string> outputTypes;					// The types which are only-output
	vector<interaction> interactions;			// The interactions
	vector<self_interaction> selfInteractions;	// The self interactions
	vector<dissociation> dissociations;			// The dissociations

	Matrix<int> interactionMat;			// A matrix of interactions
	Matrix<int> selfInteractionMat;		// A matrix of self interactions - row is the input, col is the output
	Matrix<int> dissociationMat;		// A matrix of the dissociations - row is input, col is outoput. diagonal means the self interaction exists

	// Assume output can only come from one source !!!
	vector<vector<size_t> > indexedOutputs;
	vector<vector<size_t> > indexedSelfOutputs;
	vector<vector<size_t> > indexedDissociations;

	string getOutputName(int index) {
		return ((size_t)index >= chemicalTypes.size()) ? outputTypes[index - chemicalTypes.size()] : chemicalTypes[index].name;
	}

public:
	virtual ~ChemicalNetwork() {}
	ChemicalNetwork(parsed_network &pn) {

		vector<species> &types = pn.types;
		vector<interaction_desc> &unprocessedInteractions = pn.interactions;
		vector<dissociation_desc> &unprocessedDissociations = pn.dis;

		interactionMat.resize(types.size(), types.size());
		interactionMat = -1;

		selfInteractionMat.resize(types.size(), types.size());
		selfInteractionMat = -1;

		dissociationMat.resize(types.size(), types.size());
		dissociationMat = -1;

		// First check up on all the types of species
		chemicalTypes = types;

		int pos = 0;
		for (size_t i = 0; i < types.size(); i++) {
			if (indexer.find(types[i].name) != indexer.end()) {
				cerr << "Species " << types[i].name << " appears more than once";
				exit(1);
			} else {
				indexer[types[i].name] = pos;
				pos++;
			}
		}

		indexedSelfOutputs.resize(chemicalTypes.size());
		indexedOutputs.resize(chemicalTypes.size());
		indexedDissociations.resize(chemicalTypes.size());

		// Index the interactions as well
		for (size_t i = 0; i < unprocessedInteractions.size(); i++) {
			const interaction_desc& input = unprocessedInteractions[i];

			self_interaction si;
			interaction ii;

			name_index_map::iterator first = indexer.find(input.input1);
			name_index_map::iterator second = indexer.find(input.input2);

			if (first == indexer.end()) {
				cerr << "Species " << input.input1 << " is not defined";
				exit(1);
			}
			if (second == indexer.end()) {
				cerr << "Species " << input.input2 << " is not defined";
				exit(1);
			}

			// Everything is OK - fill parameters
			if (first->second == second->second) {
				// This is a self interaction
				si.input = first->second;
			} else {
				ii.input1 = first->second;
				ii.input2 = second->second;
			}

			// Currently - no support for same interaction written twice
			if (interactionMat(first->second, second->second) != -1) { 
				cerr << "Unsupported: Interaction of " << input.input1 << " + " << input.input2 << " has appeared already" << endl;
				exit(1);				
			}

			for (size_t j = 0; j < input.outputs.size(); j++) {
				name_index_map::iterator outputType = indexer.find(input.outputs[j]);

				bool isInteracting = false;
				int outputIndex = 0;
				// The output type can be a non-interacting type
				if (outputType == indexer.end()) {
					name_index_map::iterator output = no_inter_map.find(input.outputs[j]);
					if (output == no_inter_map.end()) {
						outputIndex = chemicalTypes.size() + no_inter_map.size();
						no_inter_map[input.outputs[j]] = outputIndex;
						outputTypes.push_back(input.outputs[j]);
					} else {
						outputIndex = output->second;
					}
				} else {
					isInteracting = true;
					outputIndex = outputType->second;
				}
				if (first->second == second->second) {
					si.outputs.push_back(outputIndex);
					if (isInteracting) {
						indexedSelfOutputs[si.outputs[j]].push_back(selfInteractions.size());
						selfInteractionMat(first->second, si.outputs[j]) = selfInteractions.size();
					}
				} else {
					ii.outputs.push_back(outputIndex);
					if (isInteracting) {
						indexedOutputs[ii.outputs[j]].push_back(interactions.size());
					}
				}

			}
			if (first->second == second->second) {
				selfInteractions.push_back(si);
				interactionMat(first->second, first->second) = selfInteractions.size() - 1;
			} else {
				interactions.push_back(ii);
				interactionMat(ii.input1, ii.input2) = interactions.size() - 1;
				interactionMat(ii.input2, ii.input1) = interactions.size() - 1;
			}
		}

		// Index the Dissociations as well
		for (size_t i = 0; i < unprocessedDissociations.size(); i++) {
			const dissociation_desc& input = unprocessedDissociations[i];

			dissociation di;
			di.D = input.D;
			interaction ii;

			name_index_map::iterator first = indexer.find(input.input);			

			if (first == indexer.end()) {
				cerr << "Species " << input.input << " is not defined";
				exit(1);
			}

			// Currently - no support for same interaction written twice
			if (dissociationMat(first->second, first->second) != -1) { 
				cerr << "Unsupported: Dissociation of " << input.input << " has appeared already" << endl;
				exit(1);				
			}

			// Everything is OK - fill parameters
			di.input = first->second;		

			for (size_t j = 0; j < input.outputs.size(); j++) {
				name_index_map::iterator outputType = indexer.find(input.outputs[j]);

				bool isInteracting = false;
				int outputIndex = 0;
				// The output type can be a non-interacting type
				if (outputType == indexer.end()) {
					name_index_map::iterator output = no_inter_map.find(input.outputs[j]);
					if (output == no_inter_map.end()) {
						outputIndex = chemicalTypes.size() + no_inter_map.size();
						no_inter_map[input.outputs[j]] = outputIndex;
						outputTypes.push_back(input.outputs[j]);
					} else {
						outputIndex = output->second;
					}
				} else {
					isInteracting = true;
					outputIndex = outputType->second;
				}
				
				di.outputs.push_back(outputIndex);
				if (isInteracting) {
					indexedDissociations[di.outputs[j]].push_back(dissociations.size());
					dissociationMat(first->second, di.outputs[j]) = dissociations.size();
				}									
			}			
			dissociations.push_back(di);
			dissociationMat(first->second, first->second) = dissociations.size() - 1;			
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

		// Now, go over the dissociatios, and index the interactions containing them
		for (size_t i = 0; i < dissociations.size(); i++) {
			dissociation &di = dissociations[i];
			for (size_t j = 0; j < chemicalTypes.size(); j++) {
				if (di.input != j && dissociationMat(di.input, j) != -1) {
					di.locations.push_back(interactionMat(di.input, j));
				}
			}
		}

	}

	/**
	 * Read and parse a chemical network from an input stream.
	 * The file structure is as follows:
	 *
	 * Begin
	 * Species		Flux		W		SweepRate	Cutoff
	 * <name>		0.4		8			4		6
	 * ...
	 *
	 * Interaction
	 * <input1> <input2> => <output1> ... <outputn>
	 * <input1> <input2> => <output1> ... <outputn> 
	 * Dissociation
	 * <input> => <output1> ... <outputn> | <rate>
	 * End
	 *
	 * @param input The input stream
	 * @return a newly-created chemical network
	 */
	static parsed_network parseChemicalNetwork(istream& input) {
		vector<species> types;
		vector<interaction_desc> inters;
		vector<dissociation_desc> dis;

		enum state_enum {
			NONE,
			BEGIN, 
			SPECIES,
			INTERACTION, 
			DISSOCATION,
			END
		};

		state_enum state = NONE;
		string line;
		string str1, str2, str3, str4, str5;
		
		while (getline(input, line, '\n')) {
			istringstream is(line);
			is >> skipws;
			is >> str1;
			switch(state) {
				case NONE: 
				{														
					if (str1 == BEGIN_FILE_LABEL) {
						state = BEGIN;
					} else {
						cerr << "File should begin with 'Begin'" << endl;
						exit(1);
					}
					break;
				} 
				case BEGIN:
				{
					is >> str2 >> str3 >> str4 >> str5;
					if (str1 != SPECIES_LABEL || str2 != FLUX_LABEL || str3 != DESORPTION_LABEL || str4 != SWEEP_RATE_LABEL || str5 != CUTOFF_LABEL) {
						cerr << "File does not contain 'Species Flux Desorption SweepRate Cutoff' header" << endl;
						exit(1);
					} else {						
						state = SPECIES;
					}
					break;
				}
				case SPECIES:
				{							
					if (str1 == INTERACTION_LABEL) {
						state = INTERACTION; 
					} else if (str1 == DISSOCIATION_LABEL) {
						state = DISSOCATION;
					} else if (str1 == END_FILE_LABEL) {
						state = END;
					} else {							
						species sp;
						sp.name = str1;
						is >> sp.Flux >> sp.W >> sp.A >> sp.cutoff;
						types.push_back(sp);									
					}
					break;
				}
				case INTERACTION: 
				{
					if (str1 == END_FILE_LABEL) {						
						state = END;
					} else if (str1 == DISSOCIATION_LABEL) {
						state = DISSOCATION;
					} else {
						interaction_desc inter;
						is >> str2 >> str3;
						if (str3 != ARROW_LABEL) {
							cerr << "Interaction format: <input1> <input2> => <output1> ... <outputn>" << endl;
							exit(1);
						}
						inter.input1 = str1;
						inter.input2 = str2;
						while (!is.eof()) {
							is >> str4;
							inter.outputs.push_back(str4);
						}
						inters.push_back(inter);
					}
					break;
				}
				case DISSOCATION: 
				{
					if (str1 == INTERACTION_LABEL) {
						state = INTERACTION;
					} else if (str1 == END_FILE_LABEL) {	
						state = END;
					} else {
						dissociation_desc desc;
						desc.input = str1;
						is >> str2;												
						if (str2 != ARROW_LABEL) {
							cerr << "Dissociation format: <input1> => <output1> ... <outputn> | <rate>" << endl;
							exit(1);
						}
						is >> str3;
						while (str3 != SEPARATOR_LABEL) {							
							desc.outputs.push_back(str3);
							if (is.eof()) {
								break;
							}
							is >> str3;
						}			
						if (str3 == SEPARATOR_LABEL) {
							is >> desc.D;
						} else {
							cerr << "Dissociation format: <input1> => <output1> ... <outputn> | <rate>" << endl;
							exit(1);
						}			
						dis.push_back(desc);
					}
					break;
				} 
				case END: 
				{
					break;
				}
			}
		}		

		if (state != END) {
			cerr << "Error in file format ! " << endl;
			exit(1);
		}

		parsed_network pn;
		pn.interactions = inters;
		pn.types = types;
		pn.dis = dis;
		return pn;
	}
};

#endif
