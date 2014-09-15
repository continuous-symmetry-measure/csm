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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> //for strcmp,strlen
#include <ctype.h>  //for ispace
#include <math.h>   //for sqrt
#include "Molecule.h"
#include "babelAdapter.h"

extern "C" {
	#include "parseFunctions.h"
}

#include <openbabel/mol.h>
#include "elements.h"
#include <boost/log/trivial.hpp>

using namespace OpenBabel;

const int DIM = 3;
const double MINDOOUBLE = 1e-8;
const int LINE_BUFFER_SIZE = 1000;  /* maximal length of line of input */
const int DEPTH_ITERATIONS = 200;   /* maximal depth to descend to when checking similarity */

// ************************************************************
//       implementation
// ************************************************************

/*
 * allocates memory for the molecule structure
 */
Molecule::Molecule(int size) : _mass(size, 1.0), _valency(size, 0), _similar(size), _marked(size)
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
 * replace atom symbols with 'XX' - for unknown
 */
void Molecule::replaceSymbols()
{

	int i;

	// set new symbols 'XX'
	char *sym;

    // free old symbols
    for (i=0;i<_size;i++){
        if (_symbol[i]){
    		free(_symbol[i]);
        }
	}

    for (i=0;i<_size;i++){
    	sym = (char *)malloc(3 * sizeof(char) );
    	sym[0]='X';
    	sym[1]='X';
		sym[2]='\0';
		_symbol[i] = sym;
	}

}

/*
 * creates a molecule from an input data file
 */
Molecule* Molecule::create(FILE *in,FILE *err,bool replaceSym){

    int size,i;

    // read size
    if (fscanf(in,"%d",&size)!=1){
		BOOST_LOG_TRIVIAL(error) << "Input Error: Number of atoms not supplied";
        return NULL;
    }
    if (size == 0)
		return NULL;

    // allocate molecule
    Molecule* m = new Molecule(size);

    // read atoms
    for (i=0; i<size; i++){

        m->_symbol[i] = readString(in); // allocates space for symbol

        if(fscanf(in,"%lf%lf%lf",&(m->_pos[i][0]),&(m->_pos[i][1]),&(m->_pos[i][2]))!=3){
			BOOST_LOG_TRIVIAL(error) << "Input Error: Failed reading input for atom " << i+1;
			delete m;
            return NULL;
        }

    }

    if (replaceSym)
    	m->replaceSymbols();

    // read connectivity
    int atomNum,neighbour;
    char c;

    // allocate temporary buffer
    int *buff = (int*)malloc(size * sizeof(int));
    int pos;

    for (i=0; i<size; i++){

        // read and check atom number
        fscanf(in,"%d",&atomNum);
        if (atomNum != i /* +0 */ +1){
			BOOST_LOG_TRIVIAL(error) << "Input Error: Failed reading connectivity for atom " << i + 1;
            delete m;
            return NULL;
        }

        int valency = 0;
        pos = 0;

        // read neighbour numbers till newline
        while(1){

            // eat whitespace except newline
            c = getc(in);
            while ( isspace( (unsigned char)c ) && (c != '\n' ) )
                c = getc(in);

            // termination
            if (c == '\n')
                break;

            // otherwise put back read char
            ungetc(c,in);

            // read neighbour number
            if( (fscanf(in,"%d",&neighbour) !=1) || (neighbour > /* >= */ size) ) {
				BOOST_LOG_TRIVIAL(error) << "Input Error: Failed reading connectivity for atom " << i+1;
                delete m;
                return NULL;
            }

            // write neighbour to buffer
            buff[pos++] = neighbour  /* -0 */ -1;

            // increment valancy
            valency++;

        }
        m->_valency[i] = valency;

        // allocate memory for neighbours and copy
        m->_adjacent[i] = (int*)malloc(valency * sizeof(int));

        int j;
        for (j=0; j<valency; j++)
        	m->_adjacent[i][j] = buff[j];

    }

    free(buff);

    m->initSimilarity(DEPTH_ITERATIONS);

    return (m);
}

/*
 * creates a molecule from a PDB input data file
 */
Molecule* Molecule::createPDB(FILE *in,FILE *err,bool replaceSym){

    int size,i;

    // read size
    size = countAtomsPDB(in);
    if (size == 0)
		return NULL;

    // allocate molecule
    Molecule* m = new Molecule(size);

    // read atoms
    for (i=0; i<size; i++){

		// note: readAtom allocates space for symbol

		if (! readAtomPDB(in,&(m->_symbol[i]),m->_pos[i]) ){
			BOOST_LOG_TRIVIAL(error) << "Input Error: Failed reading input for atom " << i;
            delete m;
            return NULL;
        }

    }

    if (replaceSym)
    	m->replaceSymbols();

    // read connectivity
    int *neighbours = (int*)malloc(size * sizeof(int));
    int valency,curAtom;

    for (i=0; i<size; i++){

		if ( (! readConnectivityPDB(in,&valency,neighbours,size,&curAtom) ) ||
			(curAtom < 0) ){
			BOOST_LOG_TRIVIAL(error) << "Input Error: Failed reading connectivity element number " << i+1;
            delete m;
            return NULL;
        }

        m->_valency[curAtom] = valency;

        // allocate memory for neighbours and copy
        m->_adjacent[curAtom] = (int*)malloc(valency * sizeof(int));

        int j;
        for (j=0; j<valency; j++)
        	m->_adjacent[curAtom][j] = neighbours[j];

    }

    free(neighbours);

    m->initSimilarity(DEPTH_ITERATIONS);

    return (m);
}

/**
* Create a molecule from an OBMol
*
* @param obmol The OpenBabel Molecule
* @param replaceSym Whether to ignore atom names or not
*
* This function was moved from BabelAdapter.cpp
*/
Molecule* Molecule::createFromOBMol(OBMol &obmol, bool replaceSym, bool useMass) 
{
	int numAtoms = obmol.NumAtoms();
	int i, j;
	Molecule *mol = new Molecule(numAtoms);

	for (i = 0; i < numAtoms; i++) {
		OBAtom* atom = obmol.GetAtom(i + 1);
		if (atom->GetAtomicNum() <= ELEMENTS.size() && atom->GetAtomicNum() > 0) {
			mol->_symbol[i] = strdup(ELEMENTS[atom->GetAtomicNum() - 1].c_str());
		}
		else {
			mol->_symbol[i] = strdup(atom->GetType());
		}
		if (useMass) {
			mol->_mass[i] = atom->GetAtomicMass();
		}
		mol->_valency[i] = atom->GetValence();
		mol->_adjacent[i] = (int*)malloc(mol->_valency[i] * sizeof(int));
		j = 0;
		for (OBBondIterator itr = atom->BeginBonds(); itr != atom->EndBonds(); itr++) {
			// Check if it's the first or second atom that is the neighbour			
			if ((*itr)->GetBeginAtomIdx() == atom->GetIdx()) {
				mol->_adjacent[i][j] = (*itr)->GetEndAtomIdx() - 1;
			}
			else {
				mol->_adjacent[i][j] = (*itr)->GetBeginAtomIdx() - 1;
			}
			j++;
		}
		mol->_pos[i][0] = atom->GetX();
		mol->_pos[i][1] = atom->GetY();
		mol->_pos[i][2] = atom->GetZ();
	}

	if (replaceSym) 
		mol->replaceSymbols();

	mol->initSimilarity(DEPTH_ITERATIONS);

	return mol;
}


/*
 * Creates a new Molecule from selected atoms of the source Molecule (src)
 *
 * ! assumes selectedAtoms are present in src
 * ! This implies that the selectedAtomsSize is smaller or equal to source size
 *
 */
Molecule* Molecule::copy(int* selectedAtoms, int selectedAtomsSize, bool updateSimilarity )
{

	int i,j,old_i,len;

	// allocate molecule
    Molecule* dest = new Molecule(selectedAtomsSize);
    if (!dest){
    	return NULL; // allocation failed
    }

	//// allocate temporary buffers
	int *buff = (int*)malloc(dest->_size * sizeof(int));
	int *inversSelected = (int*)malloc(_size * sizeof(int));
	int *newGroups = (int*)malloc(_groupNum * sizeof(int));

	// init inversSelected - is the inverse of selectedAtoms (index->element,element->index)
	for ( i=0;  i< _size ;  i++ ) {
		inversSelected[i] = -1;	
	}
	for (i = 0; i < _groupNum; i++) {
		newGroups[i] = 0;
	}
	for ( i=0;  i< dest->_size ;  i++ )
		inversSelected[selectedAtoms[i]] = i;

	dest->_groupNum = 0;

	// main loop
	for (i=0; i<dest->_size; i++){

		old_i = selectedAtoms[i];

		// copy pos
		for (j=0; j<DIM; j++){
			dest->_pos[i][j] = _pos[old_i][j];
		}
		// Copy mass
		dest->_mass[i] = _mass[old_i];

		// copy similar
		if (newGroups[_similar[old_i] - 1] == 0) {
			dest->_groupNum++;
			newGroups[_similar[old_i] - 1] = dest->_groupNum;
		}
		dest->_similar[i] = newGroups[_similar[old_i] - 1];

		// copy symbol - allocate the same size as old
		len = strlen(_symbol[old_i]);
		dest->_symbol[i] = (char *)malloc((len+1) * sizeof(char) );
		strcpy(dest->_symbol[i],_symbol[old_i]);

		// update adjacency and valency
		int pos = -1;
		int valency = 0;

    	// for each item in src adjacency list
		for ( j=0;  j< _valency[old_i] ;  j++ ){

			// if item in selectedAtoms add to buffer
			pos = inversSelected[_adjacent[old_i][j]];
			if (pos != -1)
				buff[valency++] = pos;
		}

		dest->_valency[i] = valency;

		// allocate memory for adjacent and copy buffer
        dest->_adjacent[i] = (int*)malloc(valency * sizeof(int));
        for (j=0; j<valency; j++)
        	dest->_adjacent[i][j] = buff[j];

	}

	free(inversSelected);
	free(buff);
	free(newGroups);

	if (updateSimilarity)
		dest->initSimilarity(DEPTH_ITERATIONS);

	return(dest);
}

/*
 * breaks down atoms into similar groups by symbol and graph structure
 * depth -  is the desired maximal depth to take into account
 */
void Molecule::initSimilarity(int depth)
{

    int i,j,k,groupNum;

    groupNum = 1;
    setMarked(FALSE);

    // break into initial groups by symbol and valancy
    for (i=0; i<_size ; i++){

        if (_marked[i])
            continue;

        for (j=0; j<_size ; j++){
            if ( (_marked[j]) || (_valency[i] != _valency[j]) || (strcmp(_symbol[i],_symbol[j])!=0) )
                 continue;

             _similar[j] = groupNum;
             _marked[j] = true;
        }
        groupNum++;
    }

    int *group = (int*)malloc(_size * sizeof(int));     //temporary buffer
    int *subGroup = (int*)malloc(_size * sizeof(int));  //temporary buffer

    // iteratively refine the breakdown into groups
    for (k=0; k < depth ; k++){

    	setMarked(FALSE); // reset marked

    	for (i=0; i < _size; i++){

	    	if (_marked[i])
	    		continue;

	    	// mark self as done
	    	_marked[i] = true;

		    int len = 0; //group size

	    	// create the items in the group of i
	    	for(j=0; j < _size; j++)
	    		if ((j!=i) && (_similar[j] == _similar[i])){
	    			group[len++] = j;
	    		}

    		int subLen = 0; //newGroup size

	    	// for each item in the group check if it can be split or not
	    	for(j=0; j < len; j++){
	    		if (isSimilar(group[j],i))
	    			// mark similar item as done
	    			_marked[group[j]] = true;
	    		else
	    			// add to new subGroup
	    			subGroup[subLen++] = group[j];
	    	}

	    	// give subGroup members a new id
	    	bool updated = false;
	    	for(j=0; j < subLen; j++){
	    		updated = true;
	    		_similar[subGroup[j]] = groupNum;
	    	}

	    	if (updated)
	    		groupNum++;
		}
    }

    _groupNum = groupNum -1;

    free(group);
    free(subGroup);
}

/*
 * sets the general use marked array to state (true or false)
 */
void Molecule::setMarked(bool state)
{
    int i;
    for (i=0; i<_size; i++)
        _marked[i] = state;
}

/*
 * Atom a is similar to Atom b if for each neighbour of i, j has a similar neighbour
 */
int Molecule::isSimilar(int a,int b)
{
	int i, j;
	bool found = true;

	auto mark = std::vector<bool>(_size, false);

	// for each of i's neighbours
	for ( i=0;  i<_valency[a];  i++ ){

		found = false;

		for ( j=0;  j<_valency[b];  j++ ){
			if (mark[j])
				continue;

			if (_similar[_adjacent[a][i]] == _similar[_adjacent[b][j]]){
				found = true;
				mark[j] = true;
				break;
			}

		}

		if (!found)
			break;

	}

	return(found);
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
 * Creates a new Molecule from m by removing atoms who's symbol is in the
 * remove list
 */
Molecule* Molecule::stripAtoms(char** removeList, int removeListSize, int updateSimilarity)
{
	Molecule* newM;

	int i, j, count, hits;

    int *selected = (int*)malloc(_size * sizeof(int));

	// find atoms not in removeList
    count = 0;

	for ( i=0;  i< _size ;  i++ ){
		hits = 0;
		for ( j=0;  j< removeListSize ;  j++ ) {
			if( strcmp( _symbol[i], removeList[j] ) == 0 ) {
				hits++;
				break;
			}
		}
		if (hits == 0){
			selected[count] = i;
			count ++;
		}
	}

	// return a new Molecule copy
	newM = copy(selected,count,updateSimilarity);
	free(selected);
	return newM;
}

/*
 * Normalizes the position of atoms of the molecule
 * returns true if successful, false otherwise
 */
bool Molecule::normalizeMolecule(bool keepCenter = false){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	if (!keepCenter) {
		double mass_sum = 0;
		for(i=0; i< _size; i++){
			x_avg += _pos[i][0] * _mass[i];
			y_avg += _pos[i][1] * _mass[i];
			z_avg += _pos[i][2] * _mass[i];
			mass_sum += _mass[i];
		}
		x_avg /= (double)(mass_sum);
		y_avg /= (double)(mass_sum);
		z_avg /= (double)(mass_sum);
	}

	norm = 0.0;
	for(i=0; i< _size; i++){
		tmp = SQR(_pos[i][0]-x_avg) +
		      SQR(_pos[i][1]-y_avg) +
		      SQR(_pos[i][2]-z_avg);
		norm += tmp;
	}
	// normalize to 1 and not molecule size
	//norm = sqrt(norm / (double)m->_size);
	norm = sqrt(norm);


	if(norm < MINDOOUBLE)
		return false;

	for(i=0; i< _size; i++){
		_pos[i][0] = ((_pos[i][0] - x_avg) / norm);
		_pos[i][1] = ((_pos[i][1] - y_avg) / norm);
		_pos[i][2] = ((_pos[i][2] - z_avg) / norm);
	}

	_norm = norm;

	return true;
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
    m->setMarked(false);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = true;
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
    m->setMarked(false);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = true;
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
    m->setMarked(FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = true;
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
    m->setMarked(FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = true;
        }
        printf("\n");
    }
    printf("\n");

};

void Molecule::fillAtomicMasses()
{
	for (int i = 0; i < _size; i++)
		_mass[i] = getAtomicMass(_symbol[i]);
}