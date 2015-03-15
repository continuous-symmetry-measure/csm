/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>

using namespace std;

extern int main(int argc, char *argv[]); // Defined in mainRot.cpp

int RunCSM(int argc, char *argv[])
{
	return main(argc, argv);
}

int SayHello()
{
	cout << "Hello from C++" << endl;
	return 17;
}