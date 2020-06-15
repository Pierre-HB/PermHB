#pragma once

namespace factoriel {

	int const nbCityMax = 20;
	int const N = 20;

	unsigned long long factoriel(unsigned long long n);

	float random();

	void printMatrix(float mat[], int length, int maxLength);

	void decompositionFactoriel(unsigned long long n, int fact[]);

	void printArray(int array[], int length);

	void ech(int array[], int i, int j);

	void applyPermutation(int perm[], int k);

	void getPermutation(int perm[], unsigned long long n);

	float computeDistancePerm(int perm[], int nbCity, float* distanceMatrix);

	void showLong(unsigned long long n);
		

}