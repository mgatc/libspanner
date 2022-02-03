#ifndef LIBSPANNER_LW2004_H
#define LIBSPANNER_LW2004_H

#include <algorithm> // min, max
#include <cmath> // ceil
#include <unordered_set> // hashed adjacency list
#include <vector> // vertex containers

#include "constants.h"
#include "delaunay/DelaunayL2.h"
#include "../bdps/types.h"
#include "Utilities.h"



namespace spanner {

namespace lw2004 {

inline void createNewEdge(index_tPairSet &E,
                          const index_t i,
                          const index_t j)
{
    E.insert( makeNormalizedPair( i, j ) );
}

} // namespace lw2004

// alpha is set to pi/2
void LW2004(const bdps::input_t &in, bdps::output_t &out,
             number_t alpha = PI_OVER_TWO )
{
    using namespace lw2004;

    const index_t n = in.size();
    if (n > SIZE_T_MAX - 1 || n <= 1) return;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, PI_OVER_TWO ) );

    bdps::input_t P(in);
    std::vector<index_t> index;
    spatialSort<K>(P, index);

    //Step 1: Construct Delaunay triangulation
    DelaunayL2 T;

    //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
    std::vector<VertexHandle> handles(n);

    /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
      (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
    FaceHandle hint;
    for(size_t entry : index) {
        auto vh = T.insert(P[entry], hint);
        hint = vh->face();
        vh->info() = entry;
        handles[entry] = vh;
    }


    VertexHandle v_inf = T.infinite_vertex();

    //cout << "Step 1 is over...\n";
    // TriangulationPrinter tp(T);
    // tp.draw("del");
    //************* Step 2 ****************//

    std::vector<size_t> ordering;
    ordering.reserve(n);
    reverseLowDegreeOrdering(T,back_inserter(ordering));




    //************* Step 3 ****************//
    // In this step we assume alpha = pi/2 in order to minimize the degree
    index_tPairSet ePrime; // without set duplicate edges could be inserted (use the example down below)
    std::vector<bool> isProcessed(n, false);
    VertexHandle u_handle = v_inf;

    // Iterate through vertices by pi ordering
    for( size_t u : ordering ) {
        u_handle = handles.at(u);
        //assert( !T.is_infinite(u_handle) );
        isProcessed[u] = true;
        //if( printLog ) cout<<"\nu:"<<u<<" ";

        // Get neighbors of u
        VertexCirculator N = T.incident_vertices(u_handle);
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        // find a processed neighbor if it exists or we reach the start again
        while( T.is_infinite(--N) );
        VertexCirculator done(N); // set done to a vertex that is not infinite
        while( ( T.is_infinite(--N) || !isProcessed.at(N->info()) ) && N!=done );
        done = N;

        // Find and store sector boundaries, start with N
//        bool noProcessedNeighbors = !isProcessed.at( N->info() );
        std::vector<VertexHandle> sectorBoundaries{ N };
        while( --N != done ) {
            if( ( !T.is_infinite(N) && isProcessed.at( N->info() ) ) ) { // check for v_inf first or isProcessed will be out of range
                sectorBoundaries.push_back( N->handle() );
            }
        }
        assert( sectorBoundaries.size() <= 5 );


        // Now, compute the angles of the sectors, the number of cones in each sector,
        // and the actual angles
        std::vector<number_t> alphaReal( sectorBoundaries.size() );
        std::vector< std::vector<VertexHandle> > closest(sectorBoundaries.size() );

        for( size_t i=0; i<sectorBoundaries.size(); ++i ) {
            number_t sectorAngle = sectorBoundaries.size() == 1 ?
                    360.0
                    : getAngle(sectorBoundaries.at(i)->point(),
                            u_handle->point(),
                            sectorBoundaries.at((i + 1) % sectorBoundaries.size())->point() );
            auto numCones = cone_t(rint( ceil( sectorAngle / alpha ) ));
            //assert( numCones > 0 ); // guard against /0
            alphaReal[i] = sectorAngle / number_t(numCones);
            closest.at(i).resize( numCones, v_inf );
        }

        VertexHandle lastN = v_inf;
        //if( isProcessed.at( N->info() ) ) --N; // if N is processed, step
        int sector = -1;

        do { // Loop through neighbors and add appropriate edges
            if( !T.is_infinite(N) ) {
                // If N is the next sector boundary, increment sector
                if( N == sectorBoundaries.at((sector+1)%sectorBoundaries.size()) )
                    ++sector;

                if( !isProcessed[N->info()] ) {
                    //assert( sector < sectorBoundaries.size() );
                    // evaluate possible forward edges
                    number_t theta = getAngle(
                            sectorBoundaries.at(sector)->point(),
                            u_handle->point(),
                            N->point()
                    );
                    auto cone = cone_t( theta / alphaReal.at(sector) );
//                    if( cone >= closest.at(sector).size() )
//                        cone = 0;
                    // Store value until after all neighbors are processed, then add
                    if( T.is_infinite( closest.at(sector).at(cone) )
                      || CGAL::squared_distance( u_handle->point(), N->point() )
                       < CGAL::squared_distance( u_handle->point(), closest.at(sector).at(cone)->point() ) ) {
                            closest.at(sector).at(cone) = N->handle();   // if the saved vertex is infinite or longer than the current one, update
//                            if( printLog ) cout<<"s_closest["<<sector<<"]["<<cone<<"]:"<<N->info()<<" ";
                    }
                    // cross edges
                    if( !T.is_infinite( lastN ) && !isProcessed.at( lastN->info() ) ) {
//                        if( printLog ) cout<<"cross_";
                        createNewEdge(ePrime, lastN->info(), N->info());
                    }
                }
            }
            lastN = N->handle();
        } while( --N != done );

        // If N and lastN are not processed, add final cross edge
        if( !T.is_infinite(     N ) && !isProcessed.at(     N->info() )
         && !T.is_infinite( lastN ) && !isProcessed.at( lastN->info() ) )
        {
//            if( printLog ) cout<<"cross_";
            createNewEdge(ePrime, lastN->info(), N->info());
        }

        // Add edges in closest
        for( const auto& segment : closest )
            for( const auto& v : segment )
                if( !T.is_infinite(v) ) {
//                    if( printLog ) cout<<"forward_";
                    createNewEdge(ePrime, u, v->info());
                }
    }

    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    std::copy( ePrime.begin(), ePrime.end(), std::back_inserter(out) );
//    for( index_tPair e : ePrime ) {
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
//        printer.print( "lw2004" );
//        cout<<"\n";
//    }




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function LW2004

} // namespace spanner

#endif // SPANNERS_LW2004_H
