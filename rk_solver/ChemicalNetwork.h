
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

// TODO what about self interactions's output
// TODO feedback loops
// TODO more than one output types

using namespace std;

struct species {
	string name;			// The species name
	double Flux;			// Incoming Flux
	double W;		// W Rate
	double A;				// Sweeping Rate
	size_t    cutoff;			// Cutoff (max number of particles)
};

struct interaction_desc {
	string input2;			// The first species' name
	string input1;			// The second species' name
	vector<string> outputs;			// The output material
};

struct parsed_network {
	vector<species> types;
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

	typedef map<string, int> name_index_map;

	name_index_map indexer;				// A converter from a species' name to its index in the
										// chemicalTypes and meanValueHelper vectors
	name_index_map no_inter_map;		// A map of the non-interacting types


	vector<species> chemicalTypes;				// The chemical types
	vector<string> outputTypes;					// The types which are only-output
	vector<interaction> interactions;			// The interactions
	vector<self_interaction> selfInteractions;	// The self interactions
	Matrix<int> interactionMat;			// A matrix of interactions

	// Assume output can only come from one source !!!
	vector<vector<size_t> > indexedOutputs;
	vector<vector<size_t> > indexedSelfOutputs;

	string getOutputName(int index) {
		return ((size_t)index >= chemicalTypes.size()) ? outputTypes[index - chemicalTypes.size()] : chemicalTypes[index].name;
	}

	virtual void extendInteractionInfo() {}

public:
	virtual ~ChemicalNetwork() {}
	ChemicalNetwork(parsed_network &pn) {

		vector<species> &types = pn.types;
		vector<interaction_desc> &unprocessedInteractions = pn.interactions;

		interactionMat.resize(types.size(), types.size());
		interactionMat = -1;

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

			for (i = 0; i < input.outputs.size(); i++) {
				name_index_map::iterator outputType = indexer.find(input.outputs[i]);

				bool isInteracting = false;
				int outputIndex = 0;
				// The output type can be a non-interacting type
				if (outputType == indexer.end()) {
				name_index_map::iterator output = no_inter_map.find(input.outputs[i]);
					if (output == no_inter_map.end()) {
						outputIndex = chemicalTypes.size() + no_inter_map.size();
						no_inter_map[input.outputs[i]] = outputIndex;
						outputTypes.push_back(input.outputs[i]);
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
						indexedSelfOutputs[si.outputs[i]].push_back(selfInteractions.size());
					}
				} else {
					ii.outputs.push_back(outputIndex);
					if (isInteracting) {
						indexedOutputs[ii.outputs[i]].push_back(interactions.size());
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

		// Now, go over the self interactions, and index the interactions containing them
		for (size_t i = 0; i < selfInteractions.size(); i++) {
			self_interaction &si = selfInteractions[i];
			for (size_t j = 0; j < chemicalTypes.size(); j++) {
				if (si.input != j && interactionMat(si.input, j) != -1) {
					si.locations.push_back(interactionMat(si.input, j));
				}
			}
		}

		extendInteractionInfo();
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
	 * <input1> <input2> => <output1>
	 * <input1> <input2> => <output1>
	 * End
	 *
	 * @param input The input stream
	 * @return a newly-created chemical network
	 */
	static parsed_network parseChemicalNetwork(istream& input) {
		vector<species> types;
		vector<interaction_desc> inters;
		string str1, str2, str3, str4, str5;
		input >> skipws;
		input >> str1;
		if (str1 != BEGIN_FILE_LABEL) {
			cerr << "File should begin with 'Begin'" << endl;
			exit(1);
		}
		input >> str1 >> str2 >> str3 >> str4 >> str5;
		if (str1 != SPECIES_LABEL || str2 != FLUX_LABEL || str3 != DESORPTION_LABEL || str4 != SWEEP_RATE_LABEL || str5 != CUTOFF_LABEL) {
			cerr << "File does not contain 'Species Flux Desorption SweepRate Cutoff' header" << endl;
			exit(1);
		}
		input >> str1;
		while (str1 != INTERACTION_LABEL) {
			species sp;
			sp.name = str1;
			input >> sp.Flux >> sp.W >> sp.A >> sp.cutoff;
			types.push_back(sp);
			input >> str1;
			cerr << str1 << endl;
		}

		string line;
		while (getline(input, line, '\n')) {
			istringstream si(line);
			si >> skipws;
			si >> str1;
			if (str1 != END_FILE_LABEL) {
				interaction_desc inter;
				si >> str2 >> str3;
				if (str3 != ARROW_LABEL) {
					cerr << "Interaction format: <input1> <input2> => <output1>" << endl;
					exit(1);
				}
				inter.input1 = str1;
				inter.input2 = str2;
				while (!si.eof()) {
					si >> str4;
					inter.outputs.push_back(str4);
				}
				inters.push_back(inter);
				cerr << "analyzing interaction" << endl;
			} else {
				break;
			}
		}

		parsed_network pn;
		pn.interactions = inters;
		pn.types = types;
		return pn;
	}
};

#endif
