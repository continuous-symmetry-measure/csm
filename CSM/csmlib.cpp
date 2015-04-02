/*
 * CSMLIB - a library for interfacing Python code with CSM.
 *
 */

#include "csmlib.h"
#include <iostream>

using namespace std;

extern int main(int argc, char *argv[]); // Defined in mainRot.cpp

int RunCSM(const vector<string> args)
{
	vector<char *> argv;
	for (vector<string>::const_iterator it = args.begin(); it != args.end(); it++)
		argv.push_back((char *)it->c_str());

	return main((int)argv.size(), &argv[0]);
}

int SayHello()
{
	cout << "Hello from C++" << endl;
	return 17;
}