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

int AddList(const std::vector<int> numbers)
{
	int sum = 0;
	for (vector<int>::const_iterator it = numbers.begin(); it != numbers.end(); it++)
		sum += *it;
	return sum;
}

void Print(const std::string message)
{
	cout << "C++: " << message << endl;
}

std::string GetName(const std::string first, const std::string last)
{
	return last + ", " + first;
}

std::string Concat(const std::vector<string> strings)
{
	string result = "";
	for (vector<string>::const_iterator it = strings.begin(); it != strings.end(); it++)
		result += *it + " ";
	return result;
}
