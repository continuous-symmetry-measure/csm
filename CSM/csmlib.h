/*
 * The interface of the CSM library.
 * The CSM library is going to be called from Python using Cython. The library's outside interface is going
 * to be located in this file, and call internal functions as necessary.
 *
 * This file is going to change continiously, until we end up with the core C++ calculations and a big Python codebase.
 *
 * By Itay Zandbank
 */

#ifndef CSMLIB_H
#define CSMLIB_H

// Cython works much better with a C interface, use it
extern "C"
{
	// Runs the entire CSM application
	int RunCSM(int argc, char *argv[]);
	int SayHello();
}
#endif