#ifndef SPANNERS_DELAUNAYL2_H
#define SPANNERS_DELAUNAYL2_H

#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_2.h>
#include <CGAL/utils.h> // min, max

#include "tools/Utilities.h"

namespace spanner {

using namespace std;




typedef CGAL::Triangulation_vertex_base_with_info_2<index_t, K> Vb;
typedef CGAL::Triangulation_face_base_2<K>                      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                  DelaunayL2;

typedef DelaunayL2::Vertex_handle                   VertexHandle;
typedef DelaunayL2::Vertex_circulator               VertexCirculator;
typedef DelaunayL2::Face_handle                     FaceHandle;
//typedef DelaunayL2::Finite_vertices_iterator        Finite_vertices_iterator;
//typedef DelaunayL2::Finite_edges_iterator           Finite_edges_iterator;

typedef set<VertexHandle> VertexSet;
typedef unordered_set<VertexHandle> VertexHash;
template< typename V >
using VertexMap = unordered_map< VertexHandle, V >;
using AdjacencyList = VertexMap< VertexSet >;

//template< typename N >
//using EdgeInfoMap = VertexMap< VertexMap< optional<N> > >;
//typedef CGAL::Container_from_circulator<VertexCirculator> Vertex_container;


inline VertexCirculator orientCirculator(const VertexCirculator& C, const VertexHandle& v ) {
    VertexCirculator out(C),
            done(C);
    do {
        if( out->handle() == v )
            return out;
    } while( --out != done );
    return VertexCirculator();
}

class DelaunayGraph {
  public:

    /* Data */
    DelaunayL2 m_DT;
    AdjacencyList m_E;


    /* Functions */
    DelaunayGraph() = default;
    template< typename RandomAccessIterator >
    DelaunayGraph( RandomAccessIterator pointsBegin,
                   RandomAccessIterator pointsEnd )
    {
        vector<Point> P(pointsBegin, pointsEnd);
        vector<index_t> index;
        spatialSort<K>(P, index);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        for(index_t entry : index) {
            auto vh = m_DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
        }
    }

//    template< typename RandomAccessIterator >
//    void buildFromEdgeList( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
//        vector<pair<Point,Point>> edges( edgesBegin, edgesEnd );
//        std::sort( edges.begin(), edges.end() );
//
//        for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
//            VertexHandle u = m_DT.insert(e->first);
//            VertexHandle v = m_DT.insert(e->second);
//            addEdge(u, v);
//        }
//    }

    inline void addEdge(const VertexHandle& v1, const VertexHandle& v2 ) {
        addHalfEdge(v1, v2);
        addHalfEdge(v2, v1);
    }

    inline void removeEdge(const VertexHandle& v1, const VertexHandle& v2 ) {
        removeHalfEdge(v1, v2);
        removeHalfEdge(v2, v1);
    }

    inline index_t countValidNeighbors(VertexCirculator C, const VertexHash& invalid ) const {
        VertexCirculator done(C);
        int k = 0; // count N_i path length k

        do {
            if( !contains( invalid, C ) && !m_DT.is_infinite(C) )
                k++;
        } while( ++C != done );

        return k;
    }

    inline void normalizeCirculator(VertexCirculator &C, const VertexHash& invalid ) const {
        VertexCirculator done = C;
        // Position circulator so that we are guaranteed to be on the first vertex on the path N_i
        // AlgorithmFirst, loop until the circulator reaches an invalid vertex or completes a full rotation
        while( !contains( invalid, C ) && !m_DT.is_infinite(C) && ++C != done );// cout<<v_n->point()<<"\n";
        done = C;
        // Loop until the circulator reaches a valid vertex
        while( ( contains( invalid, C ) || m_DT.is_infinite(C) ) && ++C != done );// cout<<v_n->point()<<"\n";
    }

    inline index_t n() const {
        return m_DT.number_of_vertices();
    }

//    inline index_t degree() const {
//        index_t max_d = 0;
//        for( auto& v : m_E )
//            max_d = std::max( v.second.size(), max_d );
//        return max_d;
//    }

  private:

    inline void addHalfEdge(const VertexHandle& v1, const VertexHandle& v2 ) {
        auto v1IncidentVertices = m_E.find(v1);

        if(v1IncidentVertices == m_E.end() ) // v1 not found in g
            tie(v1IncidentVertices, ignore ) = m_E.insert(make_pair(v1, VertexSet() ) );

        v1IncidentVertices->second.insert(v2);
    }

    inline void removeHalfEdge(const VertexHandle& v1, const VertexHandle& v2 ) {
        auto v1IncidentVertices = m_E.find(v1);

        // if v1 in present in _E, and v2 is present in v1, remove v2 from v1
        if( contains( m_E, v1 ) && contains(v1IncidentVertices->second, v2 ) )
            v1IncidentVertices->second.erase(v2);
    }

}; // class DelaunayGraph

} // namespace spanner

#endif // SPANNERS_DELAUNAYL2_H
