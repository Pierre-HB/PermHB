#include "factoriel.h"
#include <iostream>

namespace factoriel {

    unsigned long long factoriel(unsigned long long n) {
	    if (n <= 1)return 1;
	    else return n * factoriel(n - 1);
    }

    float random() {
        return ((float)rand()) / ((float)RAND_MAX);
    }

    void printMatrix(float mat[], int length, int maxLength) {
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                std::cout << mat[i + j * maxLength] << " ";
            }
            std::cout << std::endl;
        }
    }
    void decompositionFactoriel(unsigned long long n, int fact[]) {
        unsigned long long i = 0;
        while (n != 0) {
            fact[i] = n % (i + (unsigned long long)2);
            n /= (i + (unsigned long long)2);
            i++;
        }
    }

    void printArray(int array[], int length) {
        for (int i = 0; i < length; i++) {
            std::cout << array[i] << " ";
        }
        std::cout << std::endl;
    }

    void ech(int array[], int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    void applyPermutation(int perm[], int k) {
        //a la meme effet qu'un appelle recursif avec n-1 = k
        if (k == 0)return;
        if (k % 2 == 0) {
            int temp = perm[k - 1];
            for (int i = k - 2; i >= 0; i--) {
                perm[i + 1] = perm[i];
            }
            perm[0] = temp;
            ech(perm, 0, k - 2);
            ech(perm, 1, k - 1);
        }
        else {
            ech(perm, 0, k - 1);
        }
    }

    void getPermutation(int perm[], unsigned long long n) {
        int fact[N];// = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < N; i++)fact[i] = 0;
        for (int i = 0; i < nbCityMax; i++)perm[i] = i;
        decompositionFactoriel(n, fact);
        for (int i = N - 1; i >= 0; i--) {
            //i=0 <=> n = 2
            for (int j = 0; j < fact[i]; j++) {
                applyPermutation(perm, i + 1);
                if ((i % 2) == 0) {
                    ech(perm, j, i + 1);
                }
                else {
                    ech(perm, 0, i + 1);
                }
            }
        }        
    }
    float computeDistancePerm(int perm[], int nbCity, float* distanceMatrix) {
        float distance = 0;
        for (int i = 0; i < nbCity - 1; i++) distance += distanceMatrix[perm[i] * factoriel::nbCityMax + perm[i + 1]];
        distance += distanceMatrix[perm[nbCity - 1] * factoriel::nbCityMax + perm[0]];
        return distance;
    }

    void showLongAux1(unsigned long long n) {
        if (n >= 100)std::cout << n;
        else if (n >= 10)std::cout << "0" << n;
        else std::cout << "00" << n;
    }

    void showLong(unsigned long long n) {
        if (n >= 1000) {
            showLong(n / 1000);
            std::cout << " ";
            showLongAux1(n % 1000);
        }
        else std::cout << n;
    }
}