
#ifndef LIBSPANNER_KPX2010_H
#define LIBSPANNER_KPX2010_H

#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

#include <CGAL/algorithm.h>

#include "../constants.h"
#include "../delaunay/DelaunayL2.h"
#include "../bdps/types.h"
#include "../utilities.h"



namespace spanner {

namespace kpx2010 {

bool selectEdge(index_tPairMap &E, const VertexHandle& i, const VertexHandle& j ) {
    //assert( T.is_edge( i, j ) );
    //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";

    auto existing = E.begin();
    bool inserted = false;
    tie(existing,inserted) = E.try_emplace( makeNormalizedPair( i->info(), j->info() ), false );
    if(!inserted) existing->second = true;

    return inserted;
}

} // namespace kpx2010

void KPX2010( const bdps::input_t& in, bdps::output_t& out,
              cone_t k=14 ) {
    using namespace kpx2010;

    const index_t n = in.size();
    if (n > SIZE_T_MAX - 1 || n <= 1) return;

    // ensure k >= 14
    k = std::max( k, size_t(14) );
    const number_t alpha = 2*PI / number_t(k);

    //if(printLog) cout<<"alpha:"<<alpha<<",";

    // Construct Delaunay triangulation

    std::vector<Point> P(in);
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
    index_tPairMap G_prime; // list of potential edges, value must be true for insertion to result

    // Iterate through vertices in T
    for( auto m=T.finite_vertices_begin(); m!=T.finite_vertices_end(); ++m ) {
        //if( printLog ) cout<<"\n\nm:"<<m->info()<<" ";

        // Get neighbors of m
        VertexCirculator N = T.incident_vertices(m);
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        if( T.is_infinite(N) ) --N;

        VertexCirculator done(N);
        //if(printLog) cout<<"done:"<<done->info()<<",";

        // closest vertex in each cone
        std::vector<VertexHandle> closestInCones(k, v_inf );
        // Now, let's put the neighbors found so far into a hashed set for quick lookup
        std::unordered_set<VertexHandle> selected(k);

        do { // Loop through neighbors and consider forward edges
            if( !T.is_infinite(N) ) {
                // evaluate possible forward edges
                number_t theta = getAngle(
                        done->point(),
                        m->point(),
                        N->point()
                );
                auto cone = cone_t( theta / alpha );
                //if(printLog) cout<<"N:"<<N->info()<<",theta:"<<theta<<",cone:"<<cone<<",";

                if(T.is_infinite( closestInCones.at(cone) )
                   || getDistance(m->point(), N->point()) < getDistance(m->point(), closestInCones.at(cone)->point()) )
                {   // If we made it through all that, it's the current closestInCone!
                    closestInCones[cone] = N;
                    //if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<",";
                }
            }
        } while( --N != done );

        // We've found all the closest neighbors in each cone
        // Put them in selected
        for( auto v : closestInCones )
            if( !T.is_infinite(v) )
                selected.emplace(v);

        // Now we must find every maximal sequence of empty cones
        size_t l = 0; // size of maximal empty sequence
        size_t l_local = 0; // size of current empty sequence
        size_t offset = 0; // offset in case an empty set "wraps" through the end and start of the vector
        size_t startOfSequence = 0; // start of current empty sequence
        std::unordered_set<size_t> startOfMaximalSequences(k/2);

        for( size_t i=0; i<(k+offset); ++i ) {
            //if(printLog) cout<<"i:"<<i<<",";
            if( T.is_infinite( closestInCones.at(i%k) ) ) { // empty cone
                ++l_local;          // increment
                if( l_local > l ) {  // biggest thus far, clear old starts and update l
                    startOfMaximalSequences.clear();
                    l = l_local;
                }
                if( l_local >= l )  // place the current start in the list
                    startOfMaximalSequences.emplace( startOfSequence );
                if( i+1 == k+offset ) {  // if we're about to end but on an empty sequence, keep going
                    ++offset;
                    //if(printLog) cout<<"++offset,";
                }
                //if(printLog) cout<<"l_local:"<<l_local<<",";
            } else {                    // filled cone
                //if(printLog) cout<<"filledby:"<<closestInCones.at(i%k)->info()<<",";
                l_local = 0;                 // reset l_local
                startOfSequence = (i+1) % k; // set the start of sequence to the next i
            }
        }
//        if( printLog ) {
//            cout<<"l:"<<l<<",";
//            cout<<"num_seq:"<<startOfMaximalSequences.size()<<",";
//        }
        // loop through starts of maximal sequences and add edges for them
        for( auto start : startOfMaximalSequences ) {
            number_t startAngle = number_t(start)*alpha;
//            if( printLog ) cout << "startOfMaximalSeq:"<< start<<",";
//            if( printLog ) cout << "startAngle:"<< startAngle<<",";

            while( --N != done ); // point N to reference point
            // point N to first neighbor past the empty sequence, if it exists
            while( --N != done
              && (T.is_infinite(N)
                  || getAngle(done->point(), m->point(), N->point()) < startAngle )
            );

            VertexCirculator afterSequence(N),
                              beforeSequence(N);
            while( T.is_infinite(++beforeSequence) ); // move once CCW, if it's infinite move again

//            if( printLog ) cout << "beforeSeq:"<< beforeSequence->info() <<",";
//            if( printLog ) cout << "afterSeq:"<< afterSequence->info() <<",";

            //return;
            if( l > 1 ) {
                // select the first ceil(l/2) unselected edges CCW
                auto remainingToAdd = size_t(rint(ceil(number_t(l)/2.0)));
                //if( printLog ) cout << "CCWadds:"<< remainingToAdd<<",";

                while( remainingToAdd > 0 && ++N != afterSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                        //if( printLog ) cout << "compensating:"<< N->info() <<",";
                    }
                }

                // select the first floor(l/2) unselected edges CW
                remainingToAdd = size_t(rint(floor(number_t(l)/2.0)));
                N = beforeSequence; // move N to the neighbor before the sequence
                //if( printLog ) cout << "CWadds:"<< remainingToAdd<<",";

                while( remainingToAdd > 0 && --N != beforeSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                        //if( printLog ) cout << "compensating:"<< N->info() <<",";
                    }
                }
            } else if( l == 1 ) {
                //if( printLog ) cout << "addOne,";
                VertexHandle singleSelection = v_inf;

                // consider the first CW and CCW edges (held in beforeSequence and afterSequence)
                // if one is selected already, add the other
                if( contains(selected,beforeSequence) ^ contains(selected,afterSequence) ) {
                    singleSelection = contains(selected,beforeSequence) ? afterSequence : beforeSequence;
                }
                // otherwise, add the shorter
                else if( !( contains(selected,beforeSequence) || contains(selected,afterSequence) ) ) {
                    singleSelection =
                        CGAL::squared_distance( beforeSequence->point(), m->point() )
                        < CGAL::squared_distance( afterSequence->point(), m->point() )
                            ? beforeSequence : afterSequence;
                }
                // if we found one to select, select it!
                if( !T.is_infinite(singleSelection) ) {
                    selected.emplace( singleSelection );
                    //if( printLog ) cout << "compensating:"<< singleSelection->info() <<",";
                }
            }
        }

        //bool inserted = false;
        // now add edges from each to the current vertex (u)
        for( const auto& v : selected ) {
            if( !T.is_infinite(v) ) {
                //if( printLog ) cout<<"forward_";
                selectEdge( G_prime, m, v );
            }
        }

    }

    // Send resultant graph to output iterator
    auto result = std::back_inserter(out);
    for( auto e : G_prime ) {
        if( e.second ) { // e.second holds the bool value of whether both vertices of an edge selected the edge
            *result = e.first;
            ++result;
        }
    }

} // function KPX2010

} // namespace spanner

#endif // SPANNERS_KPX2010_H
