/* Interface to the experimental C++ library that will be called from Python using Cython */

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include <string>
#include <vector>

void HelloWorld(); // Say hello
int Add(int a, int b);  // Calculate something
int AddList(const std::vector<int> numbers);
void Print(const std::string message);
std::string GetName(const std::string first, const std::string last);
std::string Concat(const std::vector<std::string> strings);

#endif