#include "salseman512.h"
#include "factoriel.h"
#include "immintrin.h" // for AVX 
#include <thread>
#include <vector>
#include <iostream>

namespace salseman512 {

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
        stream._bestDistance = _mm512_set1_ps(factoriel::nbCityMax * 2);
        stream._bestPermNb = _mm512_setzero_si512();
        stream._worstDistance = _mm512_setzero_ps();
        stream._worstPermNb = _mm512_setzero_si512();
        stream._currentPermNb = _mm512_setzero_si512();
        stream.running = true;
        return stream;
    }

    void computePermutationSalseman(int n, int NCPUT, int commonePerm[], __m512 _constantDistance, __m512 _localDistanceMatrix[], threadStream& stream) {
        //calcul une partie de (NCPUT!)*nbUnderThreads permutations
        if (n == 1) {
            __m512 _distance = _constantDistance;
            _distance = _mm512_add_ps(_distance, _localDistanceMatrix[commonePerm[0] * factoriel::nbCityMax + commonePerm[1]]);//distance += localDistanceMatrix[commonePerm[i] * nbCityMax + commonePerm[i + 1]];
            _distance = _mm512_add_ps(_distance, _localDistanceMatrix[commonePerm[NCPUT] * factoriel::nbCityMax + commonePerm[0]]);//distance += localDistanceMatrix[commonePerm[NCPUT] * nbCityMax + commonePerm[0]];

            __mmask16 mask;
            mask = _mm512_cmp_ps_mask(_distance, stream._bestDistance, _CMP_LT_OQ);//_distance < *_bestDistance
            stream._bestDistance = _mm512_mask_blend_ps(mask, stream._bestDistance, _distance);// = _distance < *_bestDistance ? _distance : *_bestDistance
            stream._bestPermNb = _mm512_mask_blend_epi8(mask, stream._bestPermNb, stream._currentPermNb);// = _distance < *_bestDistance ? *_permNb : *_bestPermNb

            mask = _mm512_cmp_ps_mask(_distance, stream._worstDistance, _CMP_GT_OQ);//_distance > *_worstDistance
            stream._worstDistance = _mm512_mask_blend_ps(mask, stream._worstDistance, _distance);// = _distance > *_worstDistance ? _distance : *_worstDistance
            stream._worstPermNb = _mm512_mask_blend_epi8(mask, stream._worstPermNb, stream._currentPermNb);// = _distance > *_worstDistance ? *_permNb : *_worstPermNb
        }
        else {
            computePermutationSalseman(n - 1, NCPUT, commonePerm, _mm512_add_ps(_constantDistance, _localDistanceMatrix[commonePerm[n-1] * factoriel::nbCityMax + commonePerm[n]]), _localDistanceMatrix, stream);
            for (int i = 0; i < n - 1; i++) {
                if ((n % 2) == 0) {
                    factoriel::ech(commonePerm, i, n - 1);
                }
                else {
                    factoriel::ech(commonePerm, 0, n - 1);
                }
                stream._currentPermNb = _mm512_add_epi32(stream._currentPermNb, _mm512_set1_epi32(1));//permNb = permNb + 1;
                computePermutationSalseman(n - 1, NCPUT, commonePerm, _mm512_add_ps(_constantDistance, _localDistanceMatrix[commonePerm[n - 1] * factoriel::nbCityMax + commonePerm[n]]), _localDistanceMatrix, stream);
            }
        }
    }

    void computePartialSalseman(unsigned long long n, float distanceMatrix[], int nbCity, int NCPUT, threadStream& stream) {
        //initialise et calcul une partie du probleme dans un Threads, chaque Therads et divise en (NCPUT!)*nbUnderThreads partie

        for (int i = 0; i < nbUnderThreads; i++)stream.nOffset[i] = (n + i) * factoriel::factoriel(NCPUT);

        int cummunePerm[factoriel::nbCityMax];
        int currentPerm[nbUnderThreads][factoriel::nbCityMax];
        for (int i = 0; i < nbCity; i++) {
            cummunePerm[i] = i;
        }

        for (int i = 0; i < nbUnderThreads; i++) factoriel::getPermutation(currentPerm[i], stream.nOffset[i]); //currentPerm[nbCity-1] = nbCity - 1 car in fixe la dernier ville, il suffit de ne calculer que les 19! autre chemin eventuel

        __m512 _partialDistance = _mm512_setzero_ps();//puisque certaine ville seront toujours placer au m�me endroit dans les diferents chemin, ont peut calculer la distance parcourue une seul fois
        for (int i = NCPUT; i < nbCity - 1; i++) {
            _partialDistance = _mm512_add_ps(_partialDistance,
                _mm512_set_ps(
                    distanceMatrix[currentPerm[0][i] * factoriel::nbCityMax + currentPerm[0][i + 1]],
                    distanceMatrix[currentPerm[1][i] * factoriel::nbCityMax + currentPerm[1][i + 1]],
                    distanceMatrix[currentPerm[2][i] * factoriel::nbCityMax + currentPerm[2][i + 1]],
                    distanceMatrix[currentPerm[3][i] * factoriel::nbCityMax + currentPerm[3][i + 1]],
                    distanceMatrix[currentPerm[4][i] * factoriel::nbCityMax + currentPerm[4][i + 1]],
                    distanceMatrix[currentPerm[5][i] * factoriel::nbCityMax + currentPerm[5][i + 1]],
                    distanceMatrix[currentPerm[6][i] * factoriel::nbCityMax + currentPerm[6][i + 1]],
                    distanceMatrix[currentPerm[7][i] * factoriel::nbCityMax + currentPerm[7][i + 1]],
                    distanceMatrix[currentPerm[8][i] * factoriel::nbCityMax + currentPerm[8][i + 1]],
                    distanceMatrix[currentPerm[9][i] * factoriel::nbCityMax + currentPerm[9][i + 1]],
                    distanceMatrix[currentPerm[10][i] * factoriel::nbCityMax + currentPerm[10][i + 1]],
                    distanceMatrix[currentPerm[11][i] * factoriel::nbCityMax + currentPerm[11][i + 1]],
                    distanceMatrix[currentPerm[12][i] * factoriel::nbCityMax + currentPerm[12][i + 1]],
                    distanceMatrix[currentPerm[13][i] * factoriel::nbCityMax + currentPerm[13][i + 1]],
                    distanceMatrix[currentPerm[14][i] * factoriel::nbCityMax + currentPerm[14][i + 1]],
                    distanceMatrix[currentPerm[15][i] * factoriel::nbCityMax + currentPerm[15][i + 1]]
                ));
        }

        // la "ville" d'indice NCPUT (la NCPUT�me ville) repr�sente l'ennemble des ville invariante durant les calculs
        __m512 _localDistanceMatrix[factoriel::nbCityMax][factoriel::nbCityMax];
        for (int i = 0; i < NCPUT; i++) {
            for (int j = 0; j < NCPUT + 1; j++) {//ATTENTION la matrice N'EST PAS symetrique
                _localDistanceMatrix[i][j] = _mm512_set_ps(
                    distanceMatrix[currentPerm[0][i] * factoriel::nbCityMax + currentPerm[0][j]],
                    distanceMatrix[currentPerm[1][i] * factoriel::nbCityMax + currentPerm[1][j]],
                    distanceMatrix[currentPerm[2][i] * factoriel::nbCityMax + currentPerm[2][j]],
                    distanceMatrix[currentPerm[3][i] * factoriel::nbCityMax + currentPerm[3][j]],
                    distanceMatrix[currentPerm[4][i] * factoriel::nbCityMax + currentPerm[4][j]],
                    distanceMatrix[currentPerm[5][i] * factoriel::nbCityMax + currentPerm[5][j]],
                    distanceMatrix[currentPerm[6][i] * factoriel::nbCityMax + currentPerm[6][j]],
                    distanceMatrix[currentPerm[7][i] * factoriel::nbCityMax + currentPerm[7][j]],
                    distanceMatrix[currentPerm[8][i] * factoriel::nbCityMax + currentPerm[8][j]],
                    distanceMatrix[currentPerm[9][i] * factoriel::nbCityMax + currentPerm[9][j]],
                    distanceMatrix[currentPerm[10][i] * factoriel::nbCityMax + currentPerm[10][j]],
                    distanceMatrix[currentPerm[11][i] * factoriel::nbCityMax + currentPerm[11][j]],
                    distanceMatrix[currentPerm[12][i] * factoriel::nbCityMax + currentPerm[12][j]],
                    distanceMatrix[currentPerm[13][i] * factoriel::nbCityMax + currentPerm[13][j]],
                    distanceMatrix[currentPerm[14][i] * factoriel::nbCityMax + currentPerm[14][j]],
                    distanceMatrix[currentPerm[15][i] * factoriel::nbCityMax + currentPerm[15][j]]
                );
            }
        }
        for (int j = 0; j < NCPUT + 1; j++) {
            _localDistanceMatrix[NCPUT][j] = _mm512_set_ps(
                distanceMatrix[currentPerm[0][nbCity - 1] * factoriel::nbCityMax + currentPerm[0][j]],
                distanceMatrix[currentPerm[1][nbCity - 1] * factoriel::nbCityMax + currentPerm[1][j]],
                distanceMatrix[currentPerm[2][nbCity - 1] * factoriel::nbCityMax + currentPerm[2][j]],
                distanceMatrix[currentPerm[3][nbCity - 1] * factoriel::nbCityMax + currentPerm[3][j]],
                distanceMatrix[currentPerm[4][nbCity - 1] * factoriel::nbCityMax + currentPerm[4][j]],
                distanceMatrix[currentPerm[5][nbCity - 1] * factoriel::nbCityMax + currentPerm[5][j]],
                distanceMatrix[currentPerm[6][nbCity - 1] * factoriel::nbCityMax + currentPerm[6][j]],
                distanceMatrix[currentPerm[7][nbCity - 1] * factoriel::nbCityMax + currentPerm[7][j]],
                distanceMatrix[currentPerm[8][nbCity - 1] * factoriel::nbCityMax + currentPerm[8][j]],
                distanceMatrix[currentPerm[9][nbCity - 1] * factoriel::nbCityMax + currentPerm[9][j]],
                distanceMatrix[currentPerm[10][nbCity - 1] * factoriel::nbCityMax + currentPerm[10][j]],
                distanceMatrix[currentPerm[11][nbCity - 1] * factoriel::nbCityMax + currentPerm[11][j]],
                distanceMatrix[currentPerm[12][nbCity - 1] * factoriel::nbCityMax + currentPerm[12][j]],
                distanceMatrix[currentPerm[13][nbCity - 1] * factoriel::nbCityMax + currentPerm[13][j]],
                distanceMatrix[currentPerm[14][nbCity - 1] * factoriel::nbCityMax + currentPerm[14][j]],
                distanceMatrix[currentPerm[15][nbCity - 1] * factoriel::nbCityMax + currentPerm[15][j]]
            );
        }

        stream._currentPermNb = _mm512_setzero_si512();

        computePermutationSalseman(NCPUT, NCPUT, cummunePerm, _partialDistance, *_localDistanceMatrix, stream);

        for (int i = 0; i < nbUnderThreads; i++) {
            stream.totalPermCompute[i] += (*(reinterpret_cast<int*>(&stream._currentPermNb) + i)) + 1;
        }
    }

    void computePartialSalsemanThread(int offset, float distanceMatrix[], int NCPUT, int nbCity, int nbThreads, threadStream* stream) {

        //initialise les variables
        (*stream).running = true;
        (*stream)._bestDistance = _mm512_set1_ps(factoriel::nbCityMax * 2);//plus grand que parcourire les diagonal autant de fois qu'il n'y a de ville
        (*stream)._bestPermNb = _mm512_setzero_si512();
        (*stream)._worstDistance = _mm512_setzero_ps();
        (*stream)._worstPermNb = _mm512_setzero_si512();
        (*stream)._currentPermNb = _mm512_setzero_si512();
        for (int i = 0; i < nbUnderThreads; i++)(*stream).bestPerm[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).worstPerm[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).nOffset[i] = 0;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).bestDistance[i] = factoriel::nbCityMax * 2;
        for (int i = 0; i < nbUnderThreads; i++)(*stream).worstDistance[i] = 0;

        for (int n = offset * nbUnderThreads; n < factoriel::factoriel(nbCity - 1) / factoriel::factoriel(NCPUT); n += nbUnderThreads * nbThreads) {
            computePartialSalseman(n, distanceMatrix, nbCity, NCPUT, (*stream));
            for (int i = 0; i < nbUnderThreads; i++) {
                if (*(reinterpret_cast<float*>(&(*stream)._bestDistance) + i) < (*stream).bestDistance[i]) {
                    (*stream).bestDistance[i] = *(reinterpret_cast<float*>(&(*stream)._bestDistance) + i);
                    (*stream).bestPerm[i] = (*stream).nOffset[nbUnderThreads - i - 1] + *(reinterpret_cast<int*>(&(*stream)._bestPermNb) + i);//les donner sont stoquer de la derniere a la premiere dans les types vectoris�
                }
                if (*(reinterpret_cast<float*>(&(*stream)._worstDistance) + i) > (*stream).worstDistance[i]) {
                    (*stream).worstDistance[i] = *(reinterpret_cast<float*>(&(*stream)._worstDistance) + i);
                    (*stream).worstPerm[i] = (*stream).nOffset[nbUnderThreads - i - 1] + *(reinterpret_cast<int*>(&(*stream)._worstPermNb) + i);//les donner sont stoquer de la derniere a la premiere dans les types vectoris�
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
            running = false;
            for (int i = 0; i < nbThreads; i++) {
                running = running || (*streams)[i].running;
                for (int j = 0; j < nbUnderThreads; j++) {
                    std::cout << " " << (unsigned long long)(((unsigned long long)(*streams)[i].nOffset[nbUnderThreads - j - 1]) + ((unsigned long long) * (reinterpret_cast<int*>(&(*streams)[i]._currentPermNb) + j))) << " |";
                    totalPerm += (*streams)[i].totalPermCompute[j];
                }
                std::cout << "| ";
            }
            std::cout << 100 * totalPerm / maxPerm << "%" << std::endl;
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
                stream.totalPermCompute[j] += streams[i].totalPermCompute[j];
            }
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
        //show the details of a stream, detail of underThreads
        std::cout << "\n\nO==========-----==========O\n\nDetailed results:\n\n";
        std::cout << "min | max\n" << std::endl;

        for (int i = 0; i < streams.size(); i++) {
            for (int j = 0; j < nbUnderThreads; j++) {
                unsigned long long bestPermNb = streams[i].bestPerm[j];
                unsigned long long worstPermNb = streams[i].worstPerm[j];
                float bestDistance = streams[i].bestDistance[j];
                float worstDistance = streams[i].worstDistance[j];

                int bestPermutation[factoriel::nbCityMax];
                factoriel::getPermutation(bestPermutation, bestPermNb);
                float bestDistanceValide = factoriel::computeDistancePerm(bestPermutation, nbCity, distanceMatrix);//0;

                int worstPermutation[factoriel::nbCityMax];
                factoriel::getPermutation(worstPermutation, worstPermNb);
                float worstDistanceValide = factoriel::computeDistancePerm(worstPermutation, nbCity, distanceMatrix);// 0;
                std::cout << "\nThread no " << i << ", underthreads no " << j << " :" << std::endl;
                std::cout << "distances  : " << bestDistance << " | " << worstDistance << std::endl;
                std::cout << "verifiees : " << bestDistanceValide << " | " << worstDistanceValide << std::endl;
                std::cout << "permNumber : " << bestPermNb << " | " << worstPermNb << std::endl;
                std::cout << "permCompute : " << streams[i].totalPermCompute[j] << std::endl;
                std::cout << "best  permutation : ";
                factoriel::printArray(bestPermutation, nbCity);
                std::cout << "worst permutation : ";
                factoriel::printArray(worstPermutation, nbCity);
            }
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
        float bestDistanceValide = factoriel::computeDistancePerm(bestPermutation, nbCity, distanceMatrix);

        int worstPermutation[factoriel::nbCityMax];
        factoriel::getPermutation(worstPermutation, worstPermNb);
        float worstDistanceValide = factoriel::computeDistancePerm(worstPermutation, nbCity, distanceMatrix);

        unsigned long long totalPermCompute = 0;
        for (int j = 0; j < nbUnderThreads; j++)totalPermCompute += stream.totalPermCompute[j];

        std::cout << "\n\nO==========-----==========O\n\nFinal results:\n\n";
        std::cout << "min | max" << std::endl;
        std::cout << "distances  : " << bestDistance << " | " << worstDistance << std::endl;
        std::cout << "verifiees : " << bestDistanceValide << " | " << worstDistanceValide << std::endl;
        std::cout << "permNumber : " << bestPermNb << " | " << worstPermNb << std::endl;
        std::cout << "permCompute : " << totalPermCompute << std::endl;
        std::cout << "best  permutation : ";
        factoriel::printArray(bestPermutation, nbCity);
        std::cout << "worst permutation : ";
        factoriel::printArray(worstPermutation, nbCity);
        std::cout << "\nComputation time : " << time << "s" << std::endl;
        std::cout << "Computation speed : ";
        factoriel::showLong((unsigned long long)(totalPermCompute / time));
        std::cout << " permutations/s" << std::endl;
    }
}