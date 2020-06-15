#include <iostream>
#include "factoriel.h"
#include "salseman32.h"
#include "salseman64.h"
#include "salseman128.h"
#include "salseman256.h"
#include "salseman512.h"
#include "enableFeatures.h"

#include "nmmintrin.h"//SEE
#include "immintrin.h" //AVX
#include "Windows.h"//pour savoir combien de Threads propose le CPU
#include "sysinfoapi.h"//pour savoir combien de Threads propose le CPU
#include <time.h>
#include <thread>
#include <string>

int nbCity = 5;
int NCPUT = 3;
int nbEnderThreads = 1;
int nbThreads = 1;
unsigned long long maxPerm = factoriel::factoriel(nbCity - 1);
int CPUThreadsMax = 0;
int CPUEnderThreadMax = 0;
int outputSize = 30;


void generateCities(float cities[]) {
    //cree une liste de ville represente par des coordonées carthésiennes
    for (int i = 0; i < nbCity; i++) {
        cities[2*i] = factoriel::random();
        cities[2*i+1] = factoriel::random();
    }
}

void generateDistanceMatrix_(float matrix[], float cities[]) {
    //matrix[y][x] = *matrix[y*length + x];
    //a cause de la représentation des matrice multidimentionel
    //matrix[x][y] = distance (x -> y)
    for (int i = 0; i < nbCity; i++) {
        for (int j = i; j < nbCity; j++) {
            float x = cities[2*i] - cities[2*j];
            float y = cities[2*i+1] - cities[2*j+1];
            matrix[i + j * factoriel::nbCityMax] = (float)sqrt(x * x + y * y);
            matrix[j + i * factoriel::nbCityMax] = matrix[i + j * factoriel::nbCityMax];
        }
    }
}


void showCities(float cities[], int size) {
    std::vector<std::vector<std::string>> pixels(size);
    for (int i = 0; i < size; i++)pixels[i].resize(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            pixels[i][j] = "  ";
        }
    }
    for (int k = 0; k < nbCity; k++) {
        int i = size * cities[2 * k];
        int j = size * cities[2 * k + 1];
        if (pixels[i][j] == "oo")pixels[i][j] = "**";
        if (pixels[i][j] == " o")pixels[i][j] = "oo";
        if (pixels[i][j] == "  ")pixels[i][j] = " o";
    }
    for (int i = 0; i < size + 1; i++)std::cout << "--";
    std::cout << "-" << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << "|";
        for (int j = 0; j < size; j++) {
            std::cout << pixels[i][j];
        }
        std::cout << " |" << std::endl;
    }
    for (int i = 0; i < size + 1; i++)std::cout << "--";
    std::cout << "-" << std::endl;

}

void showCities(float cities[], int perm[], int size) {
    std::vector<std::vector<std::string>> pixels(size);
    for (int i = 0; i < size; i++)pixels[i].resize(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            pixels[i][j] = "  ";
        }
    }
    for (int k = 0; k < nbCity; k++) {
        int i = size * cities[2 * k];
        int j = size * cities[2 * k + 1];
        if (pixels[i][j] == "oo")pixels[i][j] = "**";
        if (pixels[i][j] == " o")pixels[i][j] = "oo";
        if (pixels[i][j] == "  ")pixels[i][j] = " o";
    }
    for (int k = 0; k < nbCity; k++) {
        float x1 = cities[2 * perm[k]];
        float y1 = cities[2 * perm[k] + 1];
        float x2 = cities[2 * perm[((k + 1) % nbCity)]];
        float y2 = cities[2 * perm[((k+1)%nbCity)]+1];
        std::string pixel;
        if (y1 == y2)pixel = " -";
        else {
            float b;
            if(y1 > y2) b = (x1 - x2) / (y1 - y2);
            else b = (x2 - x1) / (y2 - y1);
            if (b < 0)pixel = " /";
            if (b >= 0)pixel = " \\";
            if (abs(b) < 0.5)pixel = "--";
            if (abs(b) >= 5)pixel = " |";
        }
        for (float a = 0; a <= 1; a += 1.0 / (float)size) {
            int i = size * (a * x1 + (1 - a) * x2);
            int j = size * (a * y1 + (1 - a) * y2);
            if(pixels[i][j] == "  ")pixels[i][j] = pixel;
        }
    }
    for (int i = 0; i < size + 1; i++)std::cout << "--";
    std::cout << "-" << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << "|";
        for (int j = 0; j < size; j++) {
            std::cout << pixels[i][j];
        }
        std::cout << " |" << std::endl;
    }
    for (int i = 0; i < size + 1; i++)std::cout << "--";
    std::cout << "-" << std::endl;
}



void init() {
    showEnableFeatures();
    std::cout << "PermutationHB vous permez de stresser votre CPU en lui faisant resoudre le problème du voyageur de commerce en exploitant le multiThreadin et les operation intrinsics";
    SYSTEM_INFO siSysInfo;
    GetSystemInfo(&siSysInfo);
    std::cout << "\n\nVotre processeur dispose de " << siSysInfo.dwNumberOfProcessors << " threads, il est inutile d'en utilise plus" << std::endl;//nombre de threads virtuels du CPU
    bool asking = true;
    while (asking) {
        std::cout << "Combien de threads voulez-vous utilise : ";
        std::cin >> nbThreads;
        asking = (nbThreads < 1);
    }
    std::cout << "\n\nIl est important de souligner que tout les processeurs ne sopporte pas les même version des operation intrinsics.";
    std::cout << "\nAinsi, il est possible que votre processeur ne puisse pas executer ce programe pour un nombre de sous-threads trop elever\n";
    asking = true;
    while (asking) {
        std::cout << "Combien de sous-threads voulez-vous utilise (1, 2, 4, 8, 16) : ";
        std::cin >> nbEnderThreads;
        asking = false;
        for (int i = 1; i <= 16; i *= 2)asking = asking || (nbEnderThreads == i);
        asking = not asking;
    }
    asking = true;
    std::cout << "\n\n";
    while (asking) {
        std::cout << "Combien de ville comporte le probleme (Le temps d'execution croit de maniere factoriel avec le nombre de ville ) (2-20) : ";
        std::cin >> nbCity;
        asking = (nbCity < 2) || (nbCity > 20);
    }
    NCPUT = nbCity /2;
    NCPUT = 10;
    if (NCPUT < 1)NCPUT = 1;
    if (NCPUT > 13)NCPUT = 13;
    maxPerm = factoriel::factoriel(nbCity - 1);
}

void showCities(float cities[], unsigned long long bestPermNb, unsigned long long worstPermNb, int size) {
    int perm[factoriel::nbCityMax];

    std::cout << "\n\nO==========-----==========O\n\nGraphiques :" << std::endl;
    showCities(cities, size);

    std::cout << "chemin le plus court : " << std::endl;
    factoriel::getPermutation(perm, bestPermNb);
    showCities(cities, perm, size);

    std::cout << "chemin le plus long : " << std::endl;
    factoriel::getPermutation(perm, worstPermNb);
    showCities(cities, perm, size);
}

int main() {
    init();
    srand(time(NULL)); // initialisation de rand
    float distanceMatrix[factoriel::nbCityMax][factoriel::nbCityMax];
    float cities[factoriel::nbCityMax][2];
    unsigned long long bestPermNb = 0;
    unsigned long long worstPermNb = 0;
    float bestDistance = 0;
    float worstDistance = 0;

    generateCities(*cities);
    generateDistanceMatrix_(*distanceMatrix, *cities);

    double duration;
    clock_t start;

    start = clock();
    if (nbEnderThreads == 1) {
        std::vector<salseman32::threadStream> streams;
        std::vector<std::thread> threads;

        salseman32::createPartialSalsemanThreads(*distanceMatrix, NCPUT, nbCity, nbThreads, streams, threads);
        salseman32::printStreams(&streams, maxPerm);

        for (int i = 0; i < nbThreads; i++) threads[i].join();

        duration = (clock() - start) / (double)CLOCKS_PER_SEC;

        salseman32::showStream(nbCity, *distanceMatrix, streams);
        salseman32::getExtremPerm(salseman32::fusionStreams(streams), bestPermNb, bestDistance, worstPermNb, worstDistance);
        showCities(*cities, bestPermNb, worstPermNb, outputSize);
        salseman32::showStreams(nbCity, *distanceMatrix, streams, duration);
    }
    if (nbEnderThreads == 2) {
        std::vector<salseman64::threadStream> streams;
        std::vector<std::thread> threads;

        salseman64::createPartialSalsemanThreads(*distanceMatrix, NCPUT, nbCity, nbThreads, streams, threads);
        salseman64::printStreams(&streams, maxPerm);

        for (int i = 0; i < nbThreads; i++) threads[i].join();

        duration = (clock() - start) / (double)CLOCKS_PER_SEC;

        salseman64::showStream(nbCity, *distanceMatrix, streams);
        salseman64::getExtremPerm(salseman64::fusionStreams(streams), bestPermNb, bestDistance, worstPermNb, worstDistance);
        showCities(*cities, bestPermNb, worstPermNb, outputSize);
        salseman64::showStreams(nbCity, *distanceMatrix, streams, duration);
    }
    if (nbEnderThreads == 4) {
        std::vector<salseman128::threadStream> streams;
        std::vector<std::thread> threads;

        salseman128::createPartialSalsemanThreads(*distanceMatrix, NCPUT, nbCity, nbThreads, streams, threads);
        salseman128::printStreams(&streams, maxPerm);

        for (int i = 0; i < nbThreads; i++) threads[i].join();

        duration = (clock() - start) / (double)CLOCKS_PER_SEC;

        salseman128::showStream(nbCity, *distanceMatrix, streams);
        salseman128::getExtremPerm(salseman128::fusionStreams(streams), bestPermNb, bestDistance, worstPermNb, worstDistance);
        showCities(*cities, bestPermNb, worstPermNb, outputSize);
        salseman128::showStreams(nbCity, *distanceMatrix, streams, duration);
    }
    if (nbEnderThreads == 8) {
        std::vector<salseman256::threadStream> streams;
        std::vector<std::thread> threads;

        salseman256::createPartialSalsemanThreads(*distanceMatrix, NCPUT, nbCity, nbThreads, streams, threads);
        salseman256::printStreams(&streams, maxPerm);

        for (int i = 0; i < nbThreads; i++) threads[i].join();

        duration = (clock() - start) / (double)CLOCKS_PER_SEC;

        salseman256::showStream(nbCity, *distanceMatrix, streams);
        salseman256::getExtremPerm(salseman256::fusionStreams(streams), bestPermNb, bestDistance, worstPermNb, worstDistance);
        showCities(*cities, bestPermNb, worstPermNb, outputSize);
        salseman256::showStreams(nbCity, *distanceMatrix, streams, duration);
    }
    if (nbEnderThreads == 16) {
        std::vector<salseman512::threadStream> streams;
        std::vector<std::thread> threads;

        salseman512::createPartialSalsemanThreads(*distanceMatrix, NCPUT, nbCity, nbThreads, streams, threads);
        salseman512::printStreams(&streams, maxPerm);

        for (int i = 0; i < nbThreads; i++) threads[i].join();

        duration = (clock() - start) / (double)CLOCKS_PER_SEC;

        salseman512::showStream(nbCity, *distanceMatrix, streams);
        salseman512::getExtremPerm(salseman512::fusionStreams(streams), bestPermNb, bestDistance, worstPermNb, worstDistance);
        showCities(*cities, bestPermNb, worstPermNb, outputSize);
        salseman512::showStreams(nbCity, *distanceMatrix, streams, duration);
    }

    std::cout << "\nappuyer sur 0 pour sortir : ";
    int temp;
    std::cin >> temp;
}
