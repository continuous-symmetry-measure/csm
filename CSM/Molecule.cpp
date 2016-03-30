/*
 * Author: shadi lahham
 *
 * A simple molecule structure with :
 *  - atom positions XYZ
 *  - atom symbols H|C|Na|Cl .. etc
 *  - adjacency matrix to represent connectivity 1|0
 *  - valency of each atom
 *  - similarity: d-similar atoms have the same symbol and same structure
 *                of children and grandchildren up to a depth of d
 *
 * Converted to C++ by Itay Zandbank
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> //for strcmp,strlen
#include <ctype.h>  //for ispace
#include <math.h>   //for sqrt
#include "Molecule.h"
#include "parseFunctions.h"
#include "elements.h"
#include "logging.h"

#include "csmlib.h"

using namespace std;

const double MINDOOUBLE = 1e-8;
const int DEPTH_ITERATIONS = 200;   /* maximal depth to descend to when checking similarity */
const int DIM = 3;

// ************************************************************
//       implementation
// ************************************************************

/*
 * allocates memory for the molecule structure
 */
Molecule::Molecule(size_t size) : _valency(size, 0), _similar(size), _mass(size, 1.0)
{
	int i;
    _size = size;

    //allocate positions
	_pos = (double **)malloc(size * sizeof(double*));
	for (i=0;i<size;i++)
		_pos[i] = (double *)malloc(DIM * sizeof(double));

    //allocate symbols
	_symbol = (char **)malloc(size * sizeof(char*));
	// individual strings allocated on reading

   
    //allocate adjacency
    _adjacent = (int **)malloc(size * sizeof(int*));
	// individual neighbours allocated on reading
};

/*
* free memory of the molecule structure
*/
Molecule::~Molecule()
{
	int i;

	// free positions
	for (i = 0; i<_size; i++){
		free(_pos[i]);
	}
	free(_pos);

	// free symbols
	for (i = 0; i<_size; i++){
		if (_symbol[i]){
			free(_symbol[i]);
		}
	}
	free(_symbol);

	// free adjacency
	for (i = 0; i<_size; i++){
		if (_adjacent[i]){
			free(_adjacent[i]);
		}
	}
	free(_adjacent);
};


/*
 * Creates a molecule from the python supplied atoms
 */
Molecule* Molecule::createFromPython(const python_molecule &molecule)
{
	const std::vector<python_atom> &atoms = molecule.atoms;
	Molecule *m = new Molecule(atoms.size());
	for (int i = 0; i < atoms.size(); i++)
	{
		m->_symbol[i] = strdup(atoms[i].symbol.c_str());
		int valency = atoms[i].adjacent.size();
		m->_valency[i] = valency;
		m->_adjacent[i] = (int*)malloc(valency * sizeof(int));

		for (int j = 0; j < valency; j++)
			m->_adjacent[i][j] = atoms[i].adjacent[j];

		for (int j = 0; j < 3; j++)
			m->_pos[i][j] = atoms[i].pos[j];

		m->_mass[i] = atoms[i].mass;
	}

	// Copy the equivalenceClasses
	if (molecule.equivalenceClasses.size())
	{
		m->_similar.resize(atoms.size(), 0);
		for (int group = 0; group < molecule.equivalenceClasses.size(); group++)
		{
			const vector<int> &eclass = molecule.equivalenceClasses[group];
			for (int i = 0; i < eclass.size(); i++)
			{
				int atom = eclass[i];
				m->_similar[atom] = group + 1;
			}
		}
		m->_groupNum = molecule.equivalenceClasses.size();
	}
	else
	{
		// Trivial equivalence classes - every atom on its own
		m->_similar.resize(atoms.size(), 0);
		for (int i = 0; i < atoms.size(); i++)
			m->_similar[i] = i;
		m->_groupNum = atoms.size(); 
	}

	return m;
}

/*
 * returns the number of elements in group 'num'
 * if there is no such group returns zero
 */
int Molecule::getGroupSize(int num)
{
	int i,j;
	i = 0;
	for (j=0; j<_size ; j++){
		if (_similar[j] == num){
			i++;
		}
	}
	return i;
}

/*
 * retrieves the similarity group number 'num' into the supplied buffer
 * returns the number of elements in that group
 * if there is no such group returns zero (buffer irrelevant)
 */
int Molecule::getGroup(int num,int* buff){
	int i,j;
	i = 0;
	for (j=0; j<_size ; j++){
		if (_similar[j] == num){
			buff[i] = j;
			i++;
		}
	}
	return i;
}

/*
 * retrieves the size of the largest similarity group
 */
int Molecule::getMaxGroupSize()
{
	int i,max=0;

	// alloc and init counting buffer
    int *count = (int*)malloc(_groupNum * sizeof(int));
    for ( i=0;  i<_groupNum;  i++ )
    	count[i] = 0;

    // count the number of members of each group
	for ( i=0;  i< _size ;  i++ ){
		count[_similar[i] - 1]++;
	}

	// get the maximum
	for ( i=0;  i< _groupNum ;  i++ ){
		max = count[i] > max ? count[i] : max;
	}

	free(count);
	return max;
}


/*
 * prints the molecule
 */
void Molecule::print()
{
    int i,j;
	Molecule *m = this; // This was previously a function accepting Molecule m as an argument
    // print molecule
    printf("molecule:\n");

    // print size
    printf("size = %d\n",m->_size);

    // print symbols and positions
	for (i=0;i<_size;i++)
        printf("%s %lf %lf %lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency & valency
    printf("\nconnectivity:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\t\tTotal = %d\n",m->_valency[i]);
    }

    printf("\nFULL similar subtree:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
        for (j=0;j<m->_size;j++){
            if (m->_similar[i] == m->_similar[j])
                printf("%d ",j+1);
        }
        printf("\n");
    }

    printf("\nsimilar subtree:\n");
	vector<bool> marked(_size, false);
   	for (i=0;i<m->_size;i++){

        if (marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( marked[j] || _similar[i] != _similar[j] )
                 continue;

             printf("%d ",j+1);
             marked[j] = true;
        }
        printf("\n");
    }

    printf("\nsimilar array:\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",i+1);
    }
    printf("\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",m->_similar[i]);
    }
    printf("\n");

};

/*
 * prints the molecule - short version
 */
void Molecule::printBasic()
{
    int i,j;
	Molecule *m = this; // This was previously a function accepting Molecule m as an argument

    // print size
    printf("%d\n",m->_size);

    // print symbols and positions
	for (i=0;i<m->_size;i++)
        printf("%2s %10lf %10lf %10lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\n");
    }

};

/*
 * prints only the similar section of the Molecule
 */
void Molecule::printSimilar()
{
    int i,j;
	Molecule *m = this; // This was previously a function accepting Molecule m as an argument

    printf("similar subtree:\n");
	vector<bool> marked(_size, false);
	for (i = 0; i<m->_size; i++){

        if (marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( marked[j] || m->_similar[i] != m->_similar[j] )
                 continue;

             printf("%d ",j+1);
             marked[j] = true;
        }
        printf("\n");
    }

    printf("\nsimilar array:\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",i+1);
    }
    printf("\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",m->_similar[i]);
    }
    printf("\n");

};

/*
 * prints the Molecule with detailed information for debugging
 */
void Molecule::printDebug()
{
    int i,j;
	Molecule *m = this; // This was previously a function accepting Molecule m as an argument

	// print similarity
	printf("Equivalent Groups:\n");
	vector<bool> marked(_size, false);
	for (i = 0; i<m->_size; i++){

        if (marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             marked[j] = true;
        }
        printf("\n");
    }

	printf("\n========DEBUG INFORMATION========\n");

    // print molecule
    printf("molecule:\n");

    // print size
    printf("size = %d\n",m->_size);

    // print symbols and positions
	for (i=0;i<m->_size;i++)
        printf("%s %lf %lf %lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency & valency
    printf("\nconnectivity:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\t\tTotal = %d\n",m->_valency[i]);
    }

    printf("\nComplete Equivalent Groups for each atom:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
        for (j=0;j<m->_size;j++){
            if (m->_similar[i] == m->_similar[j])
                printf("%d ",j+1);
        }
        printf("\n");
    }

};

/*
 * for debug purposes ... prints only the similar section of the Molecule
 */
void Molecule::printDebug2()
{
    int i,j;
	Molecule *m = this; // This was previously a function accepting Molecule m as an argument

    printf("Breakdown into groups:\n");
	vector<bool> marked(_size, false);
	for (i = 0; i<m->_size; i++){

        if (marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( marked[j] || m->_similar[i] != m->_similar[j] )
                 continue;

             printf("%d ",j+1);
             marked[j] = true;
        }
        printf("\n");
    }
    printf("\n");

};
