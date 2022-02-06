#ifndef LIBSPANNER_ASTAR_H
#define LIBSPANNER_ASTAR_H

#include <optional>
#include <unordered_map>
#include <vector>

#include <boost/heap/fibonacci_heap.hpp>

#include "../geometry.h"
#include "../types.h"
#include "../utilities.h"

namespace spanner {


    template<typename VertexContainer, typename AdjacencyList>
    std::optional<number_t> AStar(VertexContainer V, AdjacencyList G_prime, index_t start, index_t goal) {
        typedef std::pair<number_t, index_t>
                DistanceIndexPair;
        typedef boost::heap::fibonacci_heap<DistanceIndexPair, boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        Heap;
        typedef Heap::handle_type
                HeapHandle;

        index_t n = V.size();
        auto startPoint = V.at(start)->point();
        auto goalPoint = V.at(goal)->point();
        EuclideanDistanceToPoint h = {V.at(goal)->point()}; // initialize heuristic functor

        Heap open;
        std::unordered_map<index_t, HeapHandle> handleToHeap(n);
        handleToHeap[start] = open.emplace(h(startPoint), start);

        //unordered_set<size_t> closed(n);
        std::vector<index_t> parents(n);

        std::vector<number_t> g(n, INF);
        g[start] = 0;

        std::vector<number_t> f(n, INF);
        f[start] = h(startPoint);

        DistanceIndexPair current = open.top(); // initialize current vertex to start
        index_t u_index = current.second;
        auto currentPoint = startPoint;
        auto neighborPoint = currentPoint;
//    cout<<"\n    A* start:"<<startPoint;
//    cout<<",";
//    cout<<" goal:"<<V.at(goal)->point();
//    cout<<",";

        do {
            current = open.top();
            open.pop();

            u_index = current.second;
            currentPoint = V.at(u_index)->point();
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
            if (u_index == goal) return std::make_optional(g.at(goal));
//        cout<<" no goal, ";
            // loop through neighbors of current
            for (size_t neighbor : G_prime.at(u_index)) {
                neighborPoint = V.at(neighbor)->point();
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
                number_t newScore = g.at(u_index)
                                    + d(currentPoint, neighborPoint);
//            cout<<"g:"<<newScore;
//            cout<<",";
                if (newScore < g.at(neighbor)) {
                    parents[neighbor] = u_index;
                    g[neighbor] = newScore;
                    f[neighbor] = g.at(neighbor) + h(neighborPoint);
                    DistanceIndexPair q = std::make_pair(f.at(neighbor), neighbor);

                    if (contains(handleToHeap, neighbor)) {
                        HeapHandle neighborHandle = handleToHeap.at(neighbor);
                        open.update(neighborHandle, q);
                        open.update(neighborHandle);
                    } else {
                        handleToHeap[neighbor] = open.push(q);
                    }
                }
            }
        } while (!open.empty());

        return std::nullopt;
    }




/*template< typename VertexContainer, typename VertexMap, typename AdjacencyList, typename Matrix, typename H>
void AStar( const VertexContainer& V, const VertexMap& vMap, AdjacencyList& G_prime, Matrix<double>& ShortestKnownPaths, const Matrix<double>& EuclideanDistances, Matrix<H>& upperBoundHandles, size_t start, size_t goal ) {
    typedef pair<double,size_t>
        DistanceIndexPair;
    typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        Heap;
    typedef Heap::handle_type
        HeapHandle;

    size_t n = V.size();
    size_t inf = numeric_limits<size_t>::max();
    Point startPoint = V.at(start)->point();
    Point goalPoint = V.at(goal)->point();
    EuclideanDistance h = { V.at(goal)->point() }; // initialize heuristic functor

    Heap open;
    unordered_map<size_t,HeapHandle> handleToHeap(n);
    handleToHeap[start] = open.emplace( h( startPoint ), start );

    //unordered_set<size_t> closed(n);
    vector<size_t> parents(n);

    vector<double>& g = ShortestKnownPaths.at(i);

    vector<double> f( n, inf );
    f[start] = h( startPoint );

    DistanceIndexPair current = open.top(); // initialize current vertex to start
    size_t u_index = current.second;
    Point currentPoint = startPoint;
    Point neighborPoint;
//    cout<<"\n    A* start:"<<startPoint;
//    cout<<",";
//    cout<<" goal:"<<V.at(goal)->point();
//    cout<<",";

    do {
        current = open.top();
        open.pop();

        u_index = current.second;
        currentPoint = V.at(u_index)->point();
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
        if( u_index == goal ) return;
//        cout<<" no goal, ";
        double t_new = 0;
        // loop through neighbors of current
        for( size_t neighbor : G_prime.at(u_index) ) {
            neighborPoint = V.at(neighbor)->point();
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
            double newScore = g.at(u_index)
                + d( currentPoint, neighborPoint );
//            cout<<"g:"<<newScore;
//            cout<<",";
            if( newScore < g.at( neighbor ) ) {
                parents[neighbor] = u_index;
                g[neighbor] = newScore;
                f[neighbor] = g.at(neighbor) + h(neighborPoint);
                DistanceIndexPair q = make_pair( f.at(neighbor), neighbor );

                // calculate the new path's stretchFactor
                t_new = newScore / EuclideanDistances.at(i).at(u_index);
                // update t_upper in stretchFactor-Heap
                auto tValue = make_pair( t_new, make_pair(i,u_index) );
                H tHandle = upperBoundHandles.at(i).at(u_index);


                if( contains( handleToHeap, neighbor ) ) {
                    HeapHandle neighborHandle = handleToHeap.at(neighbor);
                    open.update(neighborHandle,q);
                    open.update(neighborHandle);
                } else {
                    handleToHeap[neighbor] = open.push(q);
                }
            }
        }
    } while( !open.empty() );

    return;
}*/

} // namespace spanner


#endif // LIBSPANNER_ASTAR_H


