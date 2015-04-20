/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>

using namespace std;

extern int main(int argc, char *argv[]); // Defined in mainRot.cpp

int RunCSM(csm_options options)
{
	cout << "RunCSM called with options" << endl;
	cout << "usePerm is " << options.useperm << endl;
	cout << "Output filename is " << options.outFileName;

	return 1;
}

int SayHello()
{
	cout << "Hello from C++" << endl;
	return 17;
}