#include "salseman32.h"
#include "factoriel.h"
#include <iostream>
#include <thread>
#include <vector>

namespace salseman32 {


    threadStream createThreadStream() {//type permetant un monitoring en temps reel
        threadStream stream = threadStream{};
        for (int i = 0; i < nbUnderThreads; i++) {
            stream.nOffset[i] = 0;
            stream.bestPerm[i] = 0;
            stream.worstPerm[i] = 0;
            stream.bestDistance[i] = factoriel::nbCityMax * 2;
            stream.worstDistance[i] = 0;
            stream.totalPermCompute[i] = 0;
        }
        stream._bestDistance = factoriel::nbCityMax * 2;
        stream._bestPermNb = 0;
        stream._worstDistance = 0;
        stream._worstPermNb = 0;
        stream._currentPermNb = 0;
        stream.running = true;
        return stream;
    }

    void computePermutationSalseman(int n, int NCPUT, int commonePerm[], float constantDistance, float localDistanceMatrix[], threadStream& stream) {
        //calcul une partie de (NCPUT!)*nbUnderThreads permutations
        if (n == 1) {
            float distance = constantDistance;
            distance += localDistanceMatrix[commonePerm[0] * factoriel::nbCityMax + commonePerm[1]];
            distance += localDistanceMatrix[commonePerm[NCPUT] * factoriel::nbCityMax + commonePerm[0]];
            if (distance < stream._bestDistance) {
                stream._bestDistance = distance;
                stream._bestPermNb = stream._currentPermNb;
            }
            if (distance > stream._worstDistance) {
                stream._worstDistance = distance;
                stream._worstPermNb = stream._currentPermNb;
            }
        }
        else {
            computePermutationSalseman(n - 1, NCPUT, commonePerm, constantDistance + localDistanceMatrix[commonePerm[n - 1] * factoriel::nbCityMax + commonePerm[n]], localDistanceMatrix, stream);
            for (int i = 0; i < n - 1; i++) {
                if ((n % 2) == 0) {
                    factoriel::ech(commonePerm, i, n - 1);
                }
                else {
                    factoriel::ech(commonePerm, 0, n - 1);
                }
                stream._currentPermNb+=1;
                computePermutationSalseman(n - 1, NCPUT, commonePerm, constantDistance + localDistanceMatrix[commonePerm[n -1] * factoriel::nbCityMax + commonePerm[n]], localDistanceMatrix, stream);
            }
        }
    }

    void computePartialSalseman(unsigned long long n, float distanceMatrix[], int nbCity, int NCPUT, threadStream& stream) {
        //initialise et calcul une partie du probleme dans un Threads, chaque Therads et divise en (NCPUT!)*nbUnderThreads partie

        //bien que les calcul ici pourait être plus efficasse, il n'ets pas pertinant de les optimise car ces calcul represente moins de 5%
        //du temps de calcul pour les vertion a 12 ou 13 ville (NCPUT = 6) (soit 0.3s)
        //et moins de 1% du temps pour les verssion moins trivial du problème

        for (int i = 0; i < nbUnderThreads; i++)stream.nOffset[i] = (n + i) * factoriel::factoriel(NCPUT);

        int cummunePerm[factoriel::nbCityMax];
        int currentPerm[factoriel::nbCityMax];//initialiser dans le getPermutation
        for (int i = 0; i < nbCity; i++) cummunePerm[i] = i;

        factoriel::getPermutation(currentPerm, stream.nOffset[0]); //currentPerm[nbCity-1] = nbCity - 1 car in fixe la dernier ville, il suffit de ne calculer que les 19! autre chemin eventuel
        float partialDistance = 0;//puisque certaine ville seront toujours placer au même endroit dans les diferents chemin, ont peut calculer la distance parcourue une seul fois
        for (int i = NCPUT; i < nbCity - 1; i++) {
            partialDistance += distanceMatrix[currentPerm[i] * factoriel::nbCityMax + currentPerm[i + 1]];
        }


        // la "ville" d'indice NCPUT (la NCPUTéme ville) représente l'ennemble des ville invariante durant les calculs
        float localDistanceMatrix[factoriel::nbCityMax][factoriel::nbCityMax];
        for (int i = 0; i < NCPUT; i++) {
            for (int j = 0; j < NCPUT + 1; j++) {//ATTENTION la matrice N'EST PAS symetrique
                localDistanceMatrix[i][j] = distanceMatrix[currentPerm[i] * factoriel::nbCityMax + currentPerm[j]];
            }
        }
        for (int j = 0; j < NCPUT + 1; j++) localDistanceMatrix[NCPUT][j] = distanceMatrix[currentPerm[nbCity - 1] * factoriel::nbCityMax + currentPerm[j]];

        stream._currentPermNb = 0;

        //les calculs massif seront effectuer par cette fonction, elle represente au moins 95% de la masse de calcul
        computePermutationSalseman(NCPUT, NCPUT, cummunePerm, partialDistance, *localDistanceMatrix, stream);

        stream.totalPermCompute[0] += stream._currentPermNb+1;

    }

    void computePartialSalsemanThread(int offset, float distanceMatrix[], int NCPUT, int nbCity, int nbThreads, threadStream* stream) {
        
        //initialise les variables
        (*stream).running = true;
        (*stream)._bestDistance = factoriel::nbCityMax * 2;//plus grand que parcourire les diagonal autant de fois qu'il n'y a de ville
        (*stream)._bestPermNb = 0;
        (*stream)._worstDistance = 0;
        (*stream)._worstPermNb = 0;
        (*stream)._currentPermNb = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).bestPerm[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).worstPerm[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).nOffset[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).bestDistance[i] = factoriel::nbCityMax * 2;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).worstDistance[i] = 0;

        for (int n = offset * nbUnderThreads; n < factoriel::factoriel(nbCity - 1) / factoriel::factoriel(NCPUT); n += nbUnderThreads * nbThreads) {
            computePartialSalseman(n, distanceMatrix, nbCity, NCPUT, (*stream));
            for (int i = 0; i < nbUnderThreads; i++) {
                if ((*stream)._bestDistance < (*stream).bestDistance[i]) {
                    (*stream).bestDistance[i] = (*stream)._bestDistance;
                    (*stream).bestPerm[i] = (*stream).nOffset[nbUnderThreads - i - 1] + (*stream)._bestPermNb;//les donner sont stoquer de la derniere a la premiere dans les types vectorisé
                }
                if ((*stream)._worstDistance > (*stream).worstDistance[i]) {
                    (*stream).worstDistance[i] = (*stream)._worstDistance;
                    (*stream).worstPerm[i] = (*stream).nOffset[nbUnderThreads - i - 1] + (*stream)._worstPermNb;//les donner sont stoquer de la derniere a la premiere dans les types vectorisé
                }
            }
        }
        (*stream).running = false;
    }

    void createPartialSalsemanThreads(float distanceMatrix[], int NCPUT, int nbCity, int nbThreads, std::vector<threadStream>& streams, std::vector<std::thread>& threads) {
    //genere et lance les Threads
        streams.resize(nbThreads);
        threads.resize(nbThreads);
        for (int i = 0; i < nbThreads; i++) {
            streams[i] = createThreadStream();
        }
        for (int i = 0; i < nbThreads; i++) {
            threads[i] = std::thread(computePartialSalsemanThread, i, distanceMatrix, NCPUT, nbCity, nbThreads, &streams[i]);
        }
    }

    void printStreams(std::vector<threadStream>* streams, unsigned long long const& maxPerm) {
        int nbThreads = (*streams).size();
        bool running = true;
        unsigned long long totalPerm = 0;
        
        while (running) {
            totalPerm = 0;
            for (int i = 0; i < nbThreads; i++)totalPerm += (*streams)[i].totalPermCompute[0];
            running = false;
            for (int i = 0; i < nbThreads; i++) {
                running = running || (*streams)[i].running;
                std::cout << " " << (unsigned long long)(((unsigned long long)(*streams)[i].nOffset[0]) + ((unsigned long long)(*streams)[i]._currentPermNb)) << " ||";
            }
            std::cout << " " << 100*totalPerm/maxPerm << "%" << std::endl;
        }
    }

    threadStream fusionStreams(std::vector<threadStream> const& streams) {
        int nbThreads = streams.size();
        threadStream stream = createThreadStream();
        for (int i = 0; i < nbThreads; i++) {
            for (int j = 0; j < nbUnderThreads; j++) {
                if (stream.bestDistance[j] > streams[i].bestDistance[j]) {
                    stream.bestDistance[j] = streams[i].bestDistance[j];
                    stream.bestPerm[j] = streams[i].bestPerm[j];
                }
                if (stream.worstDistance[j] < streams[i].worstDistance[j]) {
                    stream.worstDistance[j] = streams[i].worstDistance[j];
                    stream.worstPerm[j] = streams[i].worstPerm[j];
                }
            }
            stream.totalPermCompute[0] += streams[i].totalPermCompute[0];
        }
        return stream;
    }

    void getExtremPerm(threadStream const& stream, unsigned long long& bestPerm, float& bestDistance, unsigned long long& worstPerm, float& worstDistance) {
        //obtient les meilleurs et pires permutations
        bestPerm = 0;
        worstPerm = 0;
        bestDistance = factoriel::nbCityMax * 2;
        worstDistance = 0;
        for (int i = 0; i < nbUnderThreads; i++) {
            if (bestDistance > stream.bestDistance[i]) {
                bestDistance = stream.bestDistance[i];
                bestPerm = stream.bestPerm[i];
            }
            if (worstDistance < stream.worstDistance[i]) {
                worstDistance = stream.worstDistance[i];
                worstPerm = stream.worstPerm[i];
            }
        }
    }

    void showStream(int nbCity, float distanceMatrix[], std::vector<threadStream> const& streams) {
        //montre les detail d'un streams (avec les sous-Threads)
        std::cout << "\n\nO==========-----==========O\n\nDetailed results:\n\n";
        std::cout << "min | max\n" << std::endl;

        for (int i = 0; i < streams.size(); i++) {
            unsigned long long bestPermNb = streams[i].bestPerm[0];
            unsigned long long worstPermNb = streams[i].worstPerm[0];
            float bestDistance = streams[i].bestDistance[0];
            float worstDistance = streams[i].worstDistance[0];

            int bestPermutation[factoriel::nbCityMax];
            factoriel::getPermutation(bestPermutation, bestPermNb);
            float bestDistanceValide = factoriel::computeDistancePerm(bestPermutation, nbCity, distanceMatrix);//0;

            int worstPermutation[factoriel::nbCityMax];
            factoriel::getPermutation(worstPermutation, worstPermNb);
            float worstDistanceValide = factoriel::computeDistancePerm(worstPermutation, nbCity, distanceMatrix);// 0;
            std::cout << "\nThread no" << i << " :" << std::endl;
            std::cout << "distances  : " << bestDistance << " | " << worstDistance << std::endl;
            std::cout << "verifiees : " << bestDistanceValide << " | " << worstDistanceValide << std::endl;
            std::cout << "permNumber : " << bestPermNb << " | " << worstPermNb << std::endl;
            std::cout << "permCompute : " << streams[i].totalPermCompute[0] << std::endl;
            std::cout << "best  permutation : ";
            factoriel::printArray(bestPermutation, nbCity);
            std::cout << "worst permutation : ";
            factoriel::printArray(worstPermutation, nbCity);
        }
    }

    void showStreams(int nbCity, float distanceMatrix[], std::vector<threadStream> const& streams, double time) {
        //montre les information majeur d'un stream (sans les sous-Threads)
        threadStream stream = fusionStreams(streams);
        unsigned long long bestPermNb = 0;
        unsigned long long worstPermNb = 0;
        float bestDistance = 0;
        float worstDistance = 0;
        getExtremPerm(stream, bestPermNb, bestDistance, worstPermNb, worstDistance);

        int bestPermutation[factoriel::nbCityMax];
        factoriel::getPermutation(bestPermutation, bestPermNb);
        float bestDistanceValide = factoriel::computeDistancePerm(bestPermutation, nbCity, distanceMatrix);//0;

        int worstPermutation[factoriel::nbCityMax];
        factoriel::getPermutation(worstPermutation, worstPermNb);
        float worstDistanceValide = factoriel::computeDistancePerm(worstPermutation, nbCity, distanceMatrix);// 0;

        std::cout << "\n\nO==========-----==========O\n\nFinal results:\n\n";
        std::cout << "min | max" << std::endl;
        std::cout << "distances  : " << bestDistance << " | " << worstDistance << std::endl;
        std::cout << "verifiees : " << bestDistanceValide << " | " << worstDistanceValide << std::endl;
        std::cout << "permNumber : " << bestPermNb << " | " << worstPermNb << std::endl;
        std::cout << "permCompute : " << stream.totalPermCompute[0] << std::endl;
        std::cout << "best  permutation : ";
        factoriel::printArray(bestPermutation, nbCity);
        std::cout << "worst permutation : ";
        factoriel::printArray(worstPermutation, nbCity);
        std::cout << "\nComputation time : " << time << "s" << std::endl;
        std::cout << "Computation speed : ";
        factoriel::showLong((unsigned long long)(stream.totalPermCompute[0] / time));
        std::cout << " permutations/s" << std::endl;
    }
}
