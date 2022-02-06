#ifndef LIBSPANNER_FLOYDWARSHALL_H
#define LIBSPANNER_FLOYDWARSHALL_H

#include <algorithm> // swap
#include <iostream>
#include <optional>
#include <utility> // pair
#include <vector>

#include <CGAL/number_utils.h> // min
#include <CGAL/squared_distance_2.h>
#include <CGAL/utils.h>

#include "../delaunay/DelaunayL2.h"
#include "../types.h"

namespace spanner {

namespace floyd_warshall {

template< typename N >
std::optional<N> min( const std::optional<N>& ij, const std::pair< std::optional<N>,std::optional<N> >& ikj ) {
//    cout<<"min start"<<endl;

    std::optional<N> newPathLength = std::nullopt;

    if( ikj.first && ikj.second ) newPathLength = { *ikj.first + *ikj.second };
    if( ij && newPathLength ) return { CGAL::min( *ij, *newPathLength ) };
    if( ij ) return ij;

    return newPathLength;
}

} // namespace floyd_warshall

void FloydWarshall( const spanner::DelaunayGraph& G,
                    const spanner::VertexMap<size_t>& index,
                    std::vector< std::vector< std::optional<number_t> > >& distances ) {
    using namespace floyd_warshall;
    size_t N = G.size();

    // Create an NxN table to hold distances.
    std::vector< std::vector< std::optional<number_t> > > dist( N, std::vector< std::optional<number_t> >( N, std::nullopt ) );
    // container constructor should initialize optionals using default constructor, aka nullopt, aka infinity

    // Set all i==j to 0 (getDistance to self)
    for( size_t i=0; i<N; ++i )
        dist.at(i).at(i) = std::make_optional( number_t(0) );

    //assert( index.size() == N );

    // Add getDistance of each edge (u,v) in G._E to dist[u][v]
    // using indices of u and v mapped in index
    for( const auto& adjacent : G.m_E ) {
        auto u = adjacent.first; // get vertex handle
        for( const auto &v : adjacent.second )
            dist.at(index.at(u)).at(index.at(v)) = std::make_optional( CGAL::sqrt( CGAL::squared_distance( u->point(), v->point() ) ) );

    }

    // Check if going through k yields a shorter path from i to j
    for( size_t k=0; k<N; ++k )
        for( size_t i=0; i<N; ++i )
            for( size_t j=0; j<N; ++j )
                dist.at(i).at(j) = floyd_warshall::min(
                    dist.at(i).at(j),
                  { dist.at(i).at(k), dist.at(k).at(j) }
                );

    // swap the addresses for array we built with the address given in parameters
    std::swap( dist, distances );
}

} // namespace spanner


#endif // LIBSPANNER_FLOYDWARSHALL_H
