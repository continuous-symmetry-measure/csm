extern "C"
{
#include "nrutil.h"
}

#include <iostream>
#include <chrono>

#include "dvector.h"
#include "dmatrix.h"

template<typename T>
void try_vector(T& vec, int start, int size)
{
	for (int i = start; i < start + size; i++)
		vec[i] = 0;

	for (int i = start; i < start + size; i++)
		vec[i]++;

	for (int i = start + size - 2; i >= start; i--)
		vec[i + 1] += vec[i];
}

template<typename T>
void loop_vectors(T& vec, int start, int size, int count)
{
	for (int i = 0; i < count; i++)
		try_vector(vec, start, size);
}

template<typename T>
void try_matrix(T& mat, int startRow, int startCol, int size)
{
	for (int i = startRow; i < startRow + size; i++)
		for (int j = startCol; j < startCol + size; j++)
			mat[i][j] = 17;

	for (int i = startRow; i < startRow + size; i++)
		for (int j = startCol; j < startCol + size; j++)
			mat[i][j]++;

	for (int i = startRow; i < startRow + size - 1; i++)
		for (int j = startCol; j < startCol + size - 1; j++)
			mat[i][j] += mat[i + 1][j + 1];
}


void main()
{
	std::cout << "Measuring vector performances" << std::endl;

	std::cout << "C dvector" << std::endl;
	const int size = 50000000;
	const int iters = 50;
	auto dvec = dvector(1, size+1);

	auto start = std::chrono::system_clock::now();
	loop_vectors(dvec, 1, size, iters);
	auto end = std::chrono::system_clock::now();
	auto duration = end - start;
	std::cout << "dvector took " << duration.count() << std::endl << std::endl;
	free_dvector(dvec, 1, size+1);

	std::cout << "C++ dvector" << std::endl;
	csm::dvector dvec2(1, size + 1);
	start = std::chrono::system_clock::now();
	loop_vectors(dvec2, 1, size, iters);
	end = std::chrono::system_clock::now();
	duration = end - start;
	std::cout << "C++ dvector took " << duration.count() << std::endl << std::endl;

	std::cout << "C matrix" << std::endl;
	const int matSize = 4500;
	auto dmat = dmatrix(1, matSize + 1, 1, matSize + 1);
	start = std::chrono::system_clock::now();
	for (int i = 0; i < iters; i++)
		try_matrix(dmat, 1, 1, matSize);
	end = std::chrono::system_clock::now();
	duration = end - start;
	std::cout << "C dmatrix took " << duration.count() << std::endl << std::endl;
	free_dmatrix(dmat, 1, matSize + 1, 1, matSize + 1);

	std::cout << "C++ matrix" << std::endl;
	auto dmat2 = csm::dmatrix(1, matSize + 1, 1, matSize + 1);
	start = std::chrono::system_clock::now();
	for (int i = 0; i < iters; i++)
		try_matrix(dmat2, 1, 1, matSize);
	end = std::chrono::system_clock::now();
	duration = end - start;
	std::cout << "C dmatrix took " << duration.count() << std::endl << std::endl;

}
