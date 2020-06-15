#pragma once
#include "factoriel.h"
#include "nmmintrin.h"// SEE.2
#include <thread>
#include <vector>

namespace salseman128 {
    const int nbUnderThreads = 4;

    struct threadStream {
        unsigned long long nOffset[nbUnderThreads];
        unsigned long long bestPerm[nbUnderThreads];
        unsigned long long worstPerm[nbUnderThreads];
        float bestDistance[nbUnderThreads];
        float worstDistance[nbUnderThreads];
        __m128 _bestDistance;
        __m128i _bestPermNb;
        __m128 _worstDistance;
        __m128i _worstPermNb;
        __m128i _currentPermNb;
        bool running;
        unsigned long long totalPermCompute[nbUnderThreads];
    };

    void createPartialSalsemanThreads(float distanceMatrix[], int NCPUT, int nbCity, int nbThreads, std::vector<threadStream>& streams, std::vector<std::thread>& threads);

    void printStreams(std::vector<threadStream>* streams, unsigned long long const& maxPerm);

    void showStream(int nbCity, float distanceMatrix[], std::vector<threadStream> const& streams);

    void showStreams(int nbCity, float distanceMatrix[], std::vector<threadStream> const& streams, double time);

    threadStream fusionStreams(std::vector<threadStream> const& streams);

    void getExtremPerm(threadStream const& stream, unsigned long long& bestPerm, float& bestDistance, unsigned long long& worstPerm, float& worstDistance);

}