//
// Created by matt on 2/5/22.
//

#ifndef LIBSPANNER_ORDERING_H
#define LIBSPANNER_ORDERING_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <queue>

#include "include/libspanner/constants.h"
#include "include/libspanner/types.h"

namespace spanner {


    template<class K>
    void spatialSort(std::vector<typename K::Point_2> &P, std::vector<index_t> &index) {
        typedef CGAL::Spatial_sort_traits_adapter_2<K,
                typename CGAL::Pointer_property_map<typename K::Point_2>::type> SearchTraits;

        index.clear();
        index.reserve(P.size());

        std::copy(boost::counting_iterator<std::size_t>(0),
                  boost::counting_iterator<std::size_t>(P.size()),
                  std::back_inserter(index));

        CGAL::spatial_sort(index.begin(),
                           index.end(),
                           SearchTraits(CGAL::make_property_map(P)));
        //cout<<"done sorting"<<endl;
    }


    template< class DelaunayTriangulation, class VertexHash, class VertexHandle>
    inline index_t incidentChords(const DelaunayTriangulation& DT, const VertexHash& onOuterFace, const VertexHandle& v_k ) {
        typedef typename DelaunayTriangulation::Vertex_circulator VertexCirculator;
        VertexCirculator N = DT.incident_vertices(v_k),
                done(N);
        size_t c = 0;
        do {
            if( contains(onOuterFace, N->handle() ) )
                ++c;
        } while( ++N != done );
        // v_k is guaranteed to be incident to two vertices in onOuterFace, its neighbors.
        // Anything >2 is a chord
        //assert( c >= 2 );
        return (c - 2);
    }

    /**
     *  Given a Delaunay Triangulation DT and an output list out, compute the canonical out of
     *  the underlying point set.
     */
    template< class DelaunayTriangulation, typename OutputIterator >
    void canonicalOrder( const DelaunayTriangulation& DT, OutputIterator out ) {
        //Timer t(",");
        typedef typename DelaunayTriangulation::Vertex_handle VertexHandle;
        typedef typename DelaunayTriangulation::Vertex_circulator VertexCirculator;
        typedef std::unordered_set<VertexHandle> VertexHash;

        VertexHash onOuterFace, complete;
        std::queue<VertexHandle> ready;
        index_t i = DT.number_of_vertices();

        std::vector<VertexHandle> ordering(i);

        VertexCirculator v_convexHull = DT.incident_vertices(DT.infinite_vertex() ), // create a circulator of the convex hull
        done(v_convexHull );

        do { // cycle through convex hull to set vertex info
            onOuterFace.insert(v_convexHull->handle() );
            ready.push(v_convexHull->handle() );
        } while(++v_convexHull != done );

        // Reserve v_1 and v_2 so we can guarantee they are on the convex hull
        for(index_t j=0; j < 2; ++j ) {
            ordering[j] = ready.front();
            ready.pop();
        }

        VertexHandle v_k = ready.front();

        while( !ready.empty() ) {
            v_k = ready.front();
            //std::cout<<v_k->point()<<" ";
            ready.pop();
//            _algoTV.addToEventQueue( v_k, 0 );
//            _algoTV.addToEventQueue( v_k, false );

            if(incidentChords(DT, onOuterFace, v_k ) > 0 ) {
                ready.push(v_k);
                //std::cout<<"requeued";
            } else {
                //std::cout<<"processed";
                ordering[--i] = v_k;
                onOuterFace.erase(v_k);
                // add all neighbors not Complete or on outer face to ready list
                // add all neighbors not Complete to on outer face
                VertexCirculator N = DT.incident_vertices(v_k);
                done = N;
                do {
                    if( !contains( complete, N->handle() ) && !DT.is_infinite( N->handle() ) ) {
                        if( !contains(onOuterFace, N->handle() ) )
                            ready.push( N->handle() );
                        onOuterFace.insert(N->handle() );
                    }
                } while( ++N != done );
                complete.insert(v_k);
            }
            //std::cout<<"\n";
        }
        std::copy(ordering.begin(), ordering.end(), out );
    }


    template<class Graph, class OutputIterator>
    void reverseLowDegreeOrdering(const Graph &T, OutputIterator out) {

        typedef boost::heap::fibonacci_heap<index_tPair, boost::heap::compare<DegreeVertexComparator>> Heap;
        typedef Heap::handle_type HeapHandle;
        typedef typename Graph::Vertex_circulator VertexCirculator;

        const index_t n = T.number_of_vertices();

        Heap H;
        std::vector<HeapHandle> handleToHeap(n);
        //vector<size_t> piIndexedByV(n);
        std::vector<index_t> ordering(n);
        std::vector<std::unordered_set<index_t>> currentNeighbors(n);

        // Initialize the vector currentNeighbors with appropriate neighbors for every vertex
        for (auto it = T.finite_vertices_begin();
             it != T.finite_vertices_end(); ++it) {
            VertexCirculator N = T.incident_vertices(it),
                    done(N);
            do {
                if (!T.is_infinite(N))
                    currentNeighbors.at(it->info()).insert(N->info());
            } while (++N != done);

            size_t degree = currentNeighbors.at(it->info()).size();
            handleToHeap[it->info()] = H.emplace(degree, it->info());
        }

        // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
        index_t i = n - 1; // start at the last valid index

        while (!H.empty()) {
            index_tPair p = H.top();
            H.pop();
            // make sure our math is correct, e.g., degree from heap key matches neighbor container size
            //assert(p.first == currentNeighbors.at(p.second).size());
            //assert(0 <= p.first && p.first <= 5); // Lemma 1

            // Erase this vertex from incidence list of neighbors and update the neighbors' key in the heap
            for (index_t neighbor : currentNeighbors.at(p.second)) {
                currentNeighbors.at(neighbor).erase(p.second);
                HeapHandle h = handleToHeap.at(neighbor);
                index_tPair q = std::make_pair(currentNeighbors.at(neighbor).size(), neighbor);
                H.update(h, q);
                H.update(h);
            }
            currentNeighbors.at(p.second).clear();
            //piIndexedByV[p.second] = i;
            ordering[i] = p.second;
            --i;
        }
        std::copy(ordering.begin(), ordering.end(), out );

//        for (auto v : ordering) {
//            *out = v;
//            ++out;
//        }
    }

}
#endif //LIBSPANNER_ORDERING_H
