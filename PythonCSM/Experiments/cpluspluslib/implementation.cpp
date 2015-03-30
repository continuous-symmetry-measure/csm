#include <iostream>
#include "interface.h"

using namespace std;

void HelloWorld()
{
	cout << "Hello from C++" << endl;
}

int Add(int a, int b)
{
	int c = a + b;
	cout << "C++: Add(" << a << ", " << b << ") ==> " << c << endl;

	return c;
}