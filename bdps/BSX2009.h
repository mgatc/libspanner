#ifndef SPANNERS_BSX2009_H
#define SPANNERS_BSX2009_H

//#include <algorithm> // min, max
#include <cmath> // ceil
#include <vector> // vertex containers

#include "constants.h"
#include "delaunay/DelaunayLinf.h"
#include "../bdps/types.h"
#include "Utilities.h"


namespace spanner {

using namespace std;

namespace bsx2009 {

inline bool createNewEdge(//const DelaunayL2& T,
                          //const vector<VertexHandle>& handles,
                          index_tPairSet &E,
                          const index_t i, const index_t j //const index_t n, bool printLog = false
                           ) {
    //assert( std::max(i,j) < n );
    //assert( T.is_edge( handles.at(i), handles.at(j) ) );
    //if( printLog ) cout<<"add:("<<i<<","<<j<<") ";

    bool inserted = false;
    tie(ignore,inserted) = E.insert( makeNormalizedPair(i,j) );
    return inserted;
}

} // namespace bsx2009

template< typename RandomAccessIterator, typename OutputIterator >
void BSX2009( RandomAccessIterator pointsBegin,
              RandomAccessIterator pointsEnd,
              OutputIterator result,
              number_t alpha = 2*PI/3) {
    using namespace bsx2009;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, 2*PI/3 ) );
    auto numCones = cone_t( rint( ceil( 2*PI / alpha ) ) );
    //assert( numCones > 0 ); // guard against /0
    auto alphaReal =2*PI / number_t(numCones);
    //cone_t FINAL_DEGREE_BOUND = 14 + numCones;

    vector<Point> P(pointsBegin, pointsEnd);
    vector<index_t> index;
    spatialSort<K>(P, index);

    //Step 1: Construct Delaunay triangulation
    DelaunayL2 T;

    //N is the number of vertices in the delaunay triangulation.
    size_t n = P.size();
    if(n > SIZE_T_MAX - 1) return;

    //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
    vector<VertexHandle> handles(n);

    /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
      (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
    FaceHandle hint;
    for(index_t entry : index) {
        auto vh = T.insert(P[entry], hint);
        hint = vh->face();
        vh->info() = entry;
        handles[entry] = vh;
    }

    VertexHandle v_inf = T.infinite_vertex();



    //************* Step 2 ****************//
    vector<index_t> ordering;
    ordering.reserve(n);
    reverseLowDegreeOrdering(T,back_inserter(ordering));


    //************* Step 3 ****************//
    index_tPairSet ePrime;
    vector<bool> isProcessed(n, false);
    //VertexHandle u_handle = v_inf;

    // Iterate through vertices
    for( index_t u : ordering ) {
        VertexHandle u_handle = handles.at(u);
        //assert( !T.is_infinite(u_handle) );
        Point u_point = u_handle->point();
        isProcessed[u] = true;
        //if( printLog ) cout<<"\nu:"<<u<<" ";

        // Get neighbors of u
        VertexCirculator N = T.incident_vertices(u_handle );
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        if( T.is_infinite(N) ) --N;

        VertexCirculator
            closest(N),
            done(N);

        // Track degree and number of processed neighbors to ensure correctness
        size_t processedNeighbors = 0;
        size_t degree = 0;
        //bool inserted = false;

        do { // Find closest unprocessed neighbor, also count processed neighbors and neighbors in ePrime
            if( !T.is_infinite(N) ) {
                if( spanners::contains(ePrime, makeNormalizedPair(u, N->info() ) ) )
                    ++degree;

                if( isProcessed.at( N->info() ) )
                    ++processedNeighbors;
                else if( CGAL::squared_distance(u_point,N->point()) < CGAL::squared_distance(u_point,closest->point()) )
                    closest = N;
            }
        } while( --N != done );

        //if( printLog ) cout<<"degree:"<<degree<<",";
        //assert( processedNeighbors <= 5 ); // Lemma 1, proof for Lemma 3
        //assert( degree <= 15 ); // Proof for Lemma 3

        // We will add a max of numCones-1 since we are guaranteed to add the closest
        // but cannot add to the two cones touching closest.
        vector<VertexHandle> closestInCones(numCones - 1, v_inf );
        closestInCones.front() = closest; // add closest to "add" list
        while( --N != closest ); // start from the closest vertex

        // Loop through neighbors and consider forward edges
        while( --N != closest ) {
            if( !T.is_infinite(N) && !isProcessed.at(N->info()) ) {
                // evaluate possible forward edges
                number_t theta = getAngle(
                        closest->point(),
                        u_point,
                        N->point()
                );
                auto cone = cone_t( (theta-EPSILON) / alphaReal );
                // trap neighbors in forbidden cones by putting them in 0 (which is already guaranteed to be closest)
                cone = ( cone < closestInCones.size() ? cone : 0 );

                if( cone > 0 // banish the forbidden cones
                    && ( T.is_infinite( closestInCones.at(cone) )
                        || CGAL::squared_distance(u_point,N->point()) < CGAL::squared_distance(u_point,closestInCones.at(cone)->point()) ) )
                {   // If we made it through all that, it's the current closestInCone!
                    closestInCones[cone] = N->handle();
                    //if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<" ";
                }
            }
        }
        // We've found all the closest neighbors in each now,
        // now add edges from each to the current vertex (u)
        for( const auto &v : closestInCones ) {
            if( !T.is_infinite(v) ) {
                //if( printLog ) cout<<"forward_";
                createNewEdge( ePrime, u, v->info());
                //degree += size_t(inserted);
                //if( printLog ) cout<<"degree:"<<degree<<",";
            }
        }

        // Loop through neighbors again and add cross edges between
        // consecutive neighbors that are NOT processed (or infinite).
        VertexHandle lastN = N;
        do {
            --N; // Increment first, then check validity
            if( !( T.is_infinite(N) || isProcessed.at(N->info()) ) ) {
                if( !( T.is_infinite(lastN) || isProcessed.at(lastN->info()) ) ) {
                    // don't add to degree for cross edges, they are not incident on u!
                    //if( printLog ) cout<<"cross_";
                    createNewEdge( ePrime, lastN->info(), N->info() );
                }
            }
            lastN = N->handle();
        } while( N != closest );

        //assert( degree <= FINAL_DEGREE_BOUND ); // Lemma 3

    } // END OF STEP 3 LOOP



    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    std::copy( ePrime.begin(), ePrime.end(), result );
//    for( index_tPair e : ePrime ) {
//        // Edge list is only needed for printing. Remove for production.
//        //edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );
//
//        *result = e;
//        ++result;
//    }


    //
    //
    // START PRINTER NONSENSE
    //
    //

//
//    if( printLog ) {
//        GraphPrinter printer(0.007);
//        GraphPrinter::OptionsList options;
//
//        options = {
//            { "color", printer.inactiveEdgeColor },
//            { "line width", to_string(printer.inactiveEdgeWidth) }
//        };
//        printer.drawEdges( T, options );
//
//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( edgeList.begin(), edgeList.end(), options );
//
//
//        options = {
//            { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
//            { "color", make_optional( printer.backgroundColor ) }, // text color
//            { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
//            { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
//        };
//        GraphPrinter::OptionsList borderOptions = {
//            { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
//            { "color", printer.activeEdgeColor }, // additional border color
//            { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
//        };
//        printer.drawVerticesWithInfo( T, options, borderOptions );
//
//        printer.print( "bsx2009" );
//        cout<<"\n";
//    }




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function BSX2009

} // namespace spanner

#endif // SPANNERS_BSX2009_H
