#ifndef LIBSPANNER_DIJKSTRA_H
#define LIBSPANNER_DIJKSTRA_H

#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>

#include <omp.h> // <-- UNCOMMENT FOR PARALLEL, ALSO LINE 98

#include "../types.h"
#include "../utilities.h"

namespace spanner {

template<typename VertexContainer, typename AdjacencyList>
void Dijkstra(const index_t i,
              const VertexContainer &V,
              const AdjacencyList &G,
              std::vector<number_t> &ShortestPaths,
              std::vector<index_t> &Parents) {

//    typedef pair<double,size_t>
//        DistanceIndexPair;
    //typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
    typedef std::map<number_t, index_t>
            Heap;
    typedef Heap::iterator
            HeapHandle;

    const index_t n = V.size();
    auto startPoint = V.at(i);

    Heap open;
    std::unordered_map<index_t, HeapHandle> handleToHeap(n);
    handleToHeap[i] = open.emplace(0, i).first;

    ShortestPaths[i] = 0;

    auto current = open.begin(); // initialize current vertex to start
    index_t u_index = current->second;
    auto currentPoint = startPoint;
    auto neighborPoint = currentPoint;
    number_t newScore = 0;

    do {
        current = open.begin();

        u_index = current->second;
        currentPoint = V[u_index];

        // loop through neighbors of current
        for (index_t neighbor : G.at(u_index)) {
            neighborPoint = V[neighbor];

            newScore = ShortestPaths[u_index]
                       + getDistance(currentPoint, neighborPoint);

            if (newScore < ShortestPaths[neighbor]) {
                Parents[neighbor] = u_index;
                ShortestPaths[neighbor] = newScore;

                if (contains(handleToHeap, neighbor)) {
                    open.erase(handleToHeap[neighbor]);
                }
                handleToHeap[neighbor] = open.emplace(ShortestPaths[neighbor], neighbor).first;
            }
        }

        open.erase(current);
    } while (!open.empty());
}


} // spanner

#endif //LIBSPANNER_DIJKSTRA_H
