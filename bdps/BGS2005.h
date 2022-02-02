#ifndef SPANNERS_BGS2005_H
#define SPANNERS_BGS2005_H

#include <algorithm> // min, swap
#include <functional> // hash
#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "printers/GraphPrinter.h"
#include "tools/DelaunayL2.h"
#include "tools/Metrics.h"

namespace spanner {

namespace bgs2005 {

struct SplitVertexSet;
//struct SplitVertexHasher;

/* Represents a single split of a split vertex where
 * first = the vertex to be split
 * second = the split vertex's s_1
 */
struct SplitVertex {
    SplitVertexSet* V = nullptr;
    VertexHandle v;
    index_t s_1;
    VertexHandle s_1_handle;
    index_t key;

    SplitVertex() : s_1(SIZE_T_MAX), key(SIZE_T_MAX) {}
    explicit SplitVertex( const VertexHandle& v )
        : v(v), s_1(SIZE_T_MAX), key(SIZE_T_MAX) {}
    SplitVertex(const VertexHandle& v, const SplitVertex& s_1 )
        : v(v), s_1( s_1.key ), s_1_handle(s_1.v), key(SIZE_T_MAX) {}
    SplitVertex(const VertexHandle& v, const index_t& s_1 )
        : v(v), s_1( s_1 ), key(SIZE_T_MAX) {}
    SplitVertex( const SplitVertex& other ) = default;

    SplitVertex& operator=( const SplitVertex& other ) { //copy assignment
        if( this != &other ) {
            V = other.V;
            v = other.v;
            key = other.key;
            s_1 = other.s_1;
            s_1_handle = other.s_1_handle;
        }
        return *this;
    }
};

ostream& operator<<( ostream& os, const SplitVertex& v ){
    os << v.v->point() << " (" << v.s_1_handle->point() << ")";
    return os;
}

//struct SplitVertexComparator {
//    bool operator()( const SplitVertex& a, const SplitVertex& b ) {
//        return a.v < b.v || ( a.v == b.v && a.s_1 < b.s_1 );
//    }
//};

// Compare two split vertices for equality
inline bool operator==( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return lhs.v == rhs.v && lhs.s_1_handle == rhs.s_1_handle;
}

// Compare two split vertices for inequality
inline bool operator!=( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return !(lhs == rhs);
}

typedef std::set< index_t > IncidentSplitVertexContainer;
typedef std::unordered_map< index_t, IncidentSplitVertexContainer > SplitVertexEdgeMap;

struct SplitVertexSet {
    VertexMap< VertexMap< index_t > > index;

    vector< SplitVertex > V;

    index_t insert( const SplitVertex& v ) {
        index_t vectorKey;
        SplitVertex v_insert = v;
        if( contains( index, v_insert.v ) && contains( index.at(v_insert.v), v_insert.s_1_handle ) ) {
            vectorKey = index[v_insert.v][v_insert.s_1_handle];
            v_insert.key = vectorKey;
            V[vectorKey] = v_insert;
        } else {
            vectorKey = insertInContainer(v_insert);
            insertInMap(v_insert);
        }
        return vectorKey;
    }
    index_t insertInContainer(SplitVertex& v ) {
        size_t vectorKey = V.size();
        v.key = vectorKey;
        v.V = this;
        V.push_back(v);
        return vectorKey;
    }
    void insertInMap(const SplitVertex& v ) {
        index[v.v].insert_or_assign( v.s_1_handle, v.key );
    }
    SplitVertex at( const index_t i ) const {
        return V.at(i);
    }
    size_t at( const SplitVertex& v ) const {
        return index.at(v.v).at( V.at(v.s_1).v );
    }
    SplitVertex& at(const VertexHandle& v, const VertexHandle& s_1 ) {
        return V.at( index.at(v).at(s_1) );
    }
};

inline void addHalfEdge(SplitVertexEdgeMap& E, const SplitVertex& a, const SplitVertex& b ) {
    E[a.key].emplace(b.key);
}

inline size_t addVertex(SplitVertexSet& V, const VertexHandle& v, const SplitVertex& s_1 ) {
    return V.insert( SplitVertex( v, s_1 ) );
}

inline void addEdge(//const DelaunayGraph& SG,
                    SplitVertexEdgeMap& E,
                    const SplitVertex& a,
                    const SplitVertex& b ) {
    //assert( SG.m_DT.is_edge( a.v, b.v ) );

    addHalfEdge(E, b, a);
    addHalfEdge(E, a, b);
}

// Print the contents of a split vertex edge list
void print( const SplitVertexSet& V, const SplitVertexEdgeMap& E ) {
    for( const auto &v1 : E ) {
        index_t k = v1.first;
        SplitVertex u = V.at(k);
        cout<< "u: "<<u.v->point()<<" s1: "<<V.at(u.s_1).v->point()<<"\n";
        for( auto v2 : v1.second )
            cout<<"  v: "<<V.at(v2).v->point()<<" s1: "<<V.at(V.at(v2).s_1).v->point()<<"\n";
        cout<<"\n";

    }
}

// Print the contents of a split vertex list
void print( const SplitVertexSet& V ) {
    for( const auto &v : V.V ) {
        cout<<"v"<<v.key<<": "<<v;
        cout<<" s1_key: "<<v.s_1<<"\n";
    }
}

//struct SplitVertexHasher {
//    std::size_t operator()( const SplitVertex& k ) const {
//        // Compute individual hash values for first, second and third
//            // http://stackoverflow.com/a/1646913/126995
//        size_t res = 17;
//         res = res * 31 + hash< VertexHandle >()( k.v );
//         res = res * 31 + hash< VertexHandle >()( k.V->at(k.s_1).v );
//        return res;
//    }
//};



namespace spanning_graph {

inline void addFirstEdge(DelaunayGraph& G, const VertexHandle& v, const VertexCirculator& C ) {
    VertexHandle v2 = C->handle();
    G.addEdge(v, v2);
}

inline void addSecondEdge(DelaunayGraph& G, const VertexHandle& v, VertexCirculator C ) {
    while( G.m_DT.is_infinite(++C) );
    VertexHandle v2 = C->handle();
    G.addEdge(v, v2);
}

inline void addLastEdge(DelaunayGraph& G,
                        const VertexHandle& v,
                        VertexCirculator C,
                        const VertexHash& isRemoved ) {
    --C;
    VertexCirculator done(C);

    while((contains(isRemoved, C ) || G.m_DT.is_infinite(C) ) && --C != done );

    VertexHandle v2 = C->handle();

    G.addEdge(v, v2);
}

inline void removeFirstEdge(DelaunayGraph& G, VertexCirculator C ) {
    VertexHandle v1 = C->handle(),
                  v2 = (++C)->handle();
    G.removeEdge(v1, v2);
}

//inline void removeSecondEdge(DelaunayGraph& G, VertexCirculator C ) {
//    VertexHandle v1 = (++C)->handle(),
//                  v2 = (++C)->handle();
//    G.removeEdge(v1, v2);
//}

inline void removeLastEdge(DelaunayGraph& G, VertexCirculator C, const VertexHash& isRemoved ) {
    --C;

    VertexCirculator done(C);

    while((contains(isRemoved, C ) || G.m_DT.is_infinite(C) ) && --C != done );

    VertexHandle v1 = C->handle(),
                     v2 = (--C)->handle();
    G.removeEdge(v1, v2);
}

} // namespace spanning_graph

void SpanningGraph( DelaunayGraph& G ) {
    using namespace spanning_graph;

    vector< VertexHandle > canonical;

    canonicalOrder( G.m_DT ,inserter( canonical, canonical.end() ) );
    //Timer timer(",");

    VertexCirculator v_n, done;
    size_t i;

    VertexHash isRemoved(canonical.begin(), canonical.end() );

    // Add first three vertices from canonical
    for( i=0; i<3; ++i ) { // Add edges of triangle
        isRemoved.erase(canonical.at(i) );
        G.addEdge(canonical.at(i), canonical.at((i + 1) % 3));
    }
    // Add the rest of the vertices from canonical
    for( ; i<canonical.size(); ++i ) {
        isRemoved.erase(canonical.at(i) );

        v_n = G.m_DT.incident_vertices( canonical.at(i) );
        done = v_n;

        G.normalizeCirculator(v_n, isRemoved);
        done = v_n;

        index_t k = G.countValidNeighbors(v_n, isRemoved);

        if( k == 2 ) {
            // remove edge between first two vertices
            removeFirstEdge(G, v_n);
            // add edge between canonical iterator and first vertex
            addFirstEdge(G, canonical.at(i), v_n);
            // add edge between canonical iterator and second vertex
            addSecondEdge(G, canonical.at(i), v_n);

        } else if( k > 2 ) {
            // remove edge between first two vertices
            removeFirstEdge(G, v_n);
            // remove edge between last two vertices
            removeLastEdge(G, v_n, isRemoved);
            // add edge between canonical iterator and first vertex
            addFirstEdge(G, canonical.at(i), v_n);
            // add edge between canonical iterator and second vertex
            addSecondEdge(G, canonical.at(i), v_n);
            // add edge between canonical iterator and last vertex
            addLastEdge(G, canonical.at(i), v_n, isRemoved);
        }
    }

    // Test assumption
    //for( auto it=G.m_E.begin(); it!=G.m_E.end(); ++it ) { // for all v_i, 1<=i<=n
        //assert( it->second.size() <= 3 );                 // |v_i| <= 3
    //}
}



namespace transform_polygon {

inline VertexHandle find_s_1Handle(const DelaunayGraph& SG,
                                   const pair<VertexHandle, VertexMap< index_t > >& unsplit,
                                   const VertexHandle& v_n ) {
    VertexHandle v_i = unsplit.first;
    VertexCirculator N = SG.m_DT.incident_vertices(v_i); // get circulator around unsplit.first
    while( (++N)->handle() != v_n );
    while( !contains( unsplit.second, N->handle() ) ) ++N; // orient to a neighbor in unsplit.second

    return N;
}

} // namespace transform_polygon

void TransformPolygon( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& P ) {
    using namespace transform_polygon;
    P.clear();
    VertexHandle v_inf = SG.m_DT.infinite_vertex();

    // Convex hull circulator
    VertexCirculator s_1Finder = SG.m_DT.incident_vertices(v_inf ),
                      v_1Finder = SG.m_DT.incident_vertices(s_1Finder );

    v_1Finder = orientCirculator(v_1Finder, v_inf );
    ++v_1Finder; // rotate once CCW

    SplitVertex s_1(s_1Finder, 0 ); // v_1 will need an s_1 to point to

    SplitVertex v_1 = V.at(addVertex(V, v_1Finder, s_1) ),
                v_i = v_1,
                v_next = v_1;
    //cout<<v_i<<"\n";

    VertexCirculator N;

    do {
        N = SG.m_DT.incident_vertices( v_i.v );
        while( (--N)->handle() != s_1.v ); // rotate N until it points to s_1

        const VertexSet& N_SG = SG.m_E.find( v_i.v )->second; // get neighbors of v_i in SG._E
        while( !contains( N_SG, --N ) ); // rotate CW until reaching a neighbor in SG
        s_1 = v_i;
        v_next = SplitVertex( N, s_1 );

        if( v_next == v_1 ) {       // if we've looped all the way around
            SplitVertex& v_ref = V.V.at(v_1.key);
            v_ref.s_1 = s_1.key; // update v_1's s_1
            v_i = v_ref; // update v_i to explicitly point to v_1
        } else { // add new vertex
            v_i = V.at(addVertex(V, N, s_1) );
        }

        //assert( !SG.m_DT.is_infinite(N) );
        addEdge( P, v_i, s_1);

    } while( v_i != v_1 );

    for( auto e=SG.m_DT.finite_edges_begin(); e!=SG.m_DT.finite_edges_end(); ++e ) {
        VertexHandle u = e->first->vertex((e->second + 1) % 3 ),
                         v = e->first->vertex( (e->second+2)%3 );
        if( !contains( SG.m_E.at(u), v ) ) { // only add edges that are not in SG
            pair< const VertexHandle, VertexMap< index_t > >&
                u_unsplit = *V.index.find(u ),
                v_unsplit = *V.index.find(v );
            VertexHandle s_1_u = find_s_1Handle(SG, u_unsplit, v),
                             s_1_v = find_s_1Handle(SG, v_unsplit, u);
            addEdge( P, V.at(u, s_1_u), V.at(v, s_1_v));
        }
    }
}



namespace polygon_spanner {

enum VertexStatus { Known, Complete };
using VertexStatusMap = unordered_map< index_t, VertexStatus >;



struct VertexHandleHash {
    size_t operator()( const SplitVertex& k ) const {
        return hash< VertexHandle >()(k.v );
    }
};
struct VertexHandleComparator {
    bool operator()( const SplitVertex& a, const SplitVertex& b ) const {
        return a.v == b.v;
    }
};

/*
 * Iterates from s_1 to s_m (inclusive) and performs "foreach" on each vertex.
 * Provides a SplitVertex of the currently pointed-to vertex.
 */
template< typename F >
void forEachNeighbor(const DelaunayGraph& SG,
                     const SplitVertexSet& V,
                     const SplitVertexEdgeMap& E,
                     const SplitVertex& v_i,
                     F foreach ) {
    //assert( !SG.m_DT.is_infinite( v_i.v ) );
    /*
     *    The circulator provided by SG._DT will provide the correct ordering.
     * However, the vertices provided by the circulator are not split. Aside
     * from s_1 and s_m, which can have equal primary vertices (with differing
     * s_1s of course), the primary vertices of these neighbors will be unique.
     * Therefore, we will create an unordered set of incident vertices to v_i
     * and provide a custom hash function that will simply pass along the
     * CGAL VertexHandle hash. This way, we can match a split vertex from the
     * vertex handle provided by the circulator in constant time.
     */
    unordered_set< SplitVertex, VertexHandleHash, VertexHandleComparator > N_E;

    for( auto k: E.at( v_i.key ) ) {
        if( k != v_i.s_1 ) { // don't enter s_1
            N_E.insert( V.at(k) );
        }
    }
    //cout<<"N_E size:"<<N_E.size()<<"\n";

    VertexCirculator N = SG.m_DT.incident_vertices(v_i.v );
    N = orientCirculator( N, v_i.s_1_handle );


    SplitVertex v_neighbor = V.at( v_i.s_1 );
    bool include_s_1 = contains( E.at(v_i.key), v_neighbor.key );
    //cout<< "    start foreach: "<<v_i<<" n: "<<(N_E.size()+size_t(include_s_1))<<"\n";

    if( include_s_1 ) {
    //    cout<<"    v_neighbor: "<<v_neighbor<<"\n";
        foreach( v_neighbor );
    }

    do {
        --N;
        v_neighbor = SplitVertex( N->handle() );
        if( contains( N_E, v_neighbor ) ) {
            v_neighbor = *N_E.find( v_neighbor );
            //cout<<"    v_neighbor: "<<v_neighbor<<"\n";
            foreach( v_neighbor );
            N_E.erase( v_neighbor );
        }
    } while( !N_E.empty() );
}

inline SplitVertex get_s_m( const DelaunayGraph& SG, const SplitVertexSet& V, const SplitVertexEdgeMap& E, const SplitVertex& v_split ) {
    index_t s_i;
    forEachNeighbor(SG, V, E, v_split, [&](SplitVertex &v_n) {
        s_i = v_n.key;
    });
    return V.at(s_i);
}

/*
 * Adds a split vertex edge, asserting that neither vertices a nor b are marked as Complete
 */
inline void addPolygonSpannerEdge(//const DelaunayGraph& SG,
                                  SplitVertexEdgeMap& E,
                                  SplitVertex& a,
                                  SplitVertex& b ) {
    //assert( !contains( status, a.key ) || status.at(a.key) != Complete );
    //assert( !contains( status, b.key ) || status.at(b.key) != Complete );
    return addEdge( E, a, b);
}

void addCrossEdges(const DelaunayGraph& SG,
                   SplitVertexSet& V,
                   SplitVertexEdgeMap& E,
                   SplitVertexEdgeMap& E_P,
                   SplitVertex& p,
                   SplitVertex& q,
                   SplitVertex& r) {
    optional< SplitVertex > v_last = nullopt;
    bool isInZone = false;
    forEachNeighbor(SG, V, E, q, [&](SplitVertex &v_n) {
        if (isInZone && v_last) {
            addPolygonSpannerEdge(E_P, *v_last, v_n);
            //SG.addToEventQueue( {s_last.first, N}, true );
        }
        if (v_n.key == p.key) isInZone = true;
        if (v_n.key == r.key) isInZone = false;

        v_last = {v_n};
    });
}

void addForwardEdges(const DelaunayGraph& SG,
                     SplitVertexSet& V,
                     SplitVertexEdgeMap& E,
                     SplitVertexEdgeMap& E_P,
                     SplitVertex& p,
                     SplitVertex& q,
                     SplitVertex& r) {

    //double deg = 180/PI; // used for displaying angles in degree

    number_t alpha = getAngle(p.v->point(), q.v->point(), r.v->point());
    auto subangles = cone_t(rint( ceil( alpha / (PI/2) ) ));
    auto beta = alpha / number_t(subangles);

    vector< SplitVertex > add( subangles, SplitVertex( SG.m_DT.infinite_vertex() ) ); // initialize add to infinite vertex

    number_t theta;
    cone_t i;
    bool isInZone = false;

    forEachNeighbor(SG, V, E, q, [&](SplitVertex &v_n) {
        // r is guaranteed to be already added or will be added immediately after this step
        // so set isInZone to false before processing it
        if (v_n.key == r.key) isInZone = false;
        //SG.addToEventQueue( N, 1 );// focus1 on N
        if (isInZone) {
            theta = getAngle(p.v->point(), q.v->point(), v_n.v->point());
            if (theta > 2 * PI - EPSILON)
                theta = 0;
            i = CGAL::min(cone_t(floor(theta / beta)), subangles - 1);

            if (SG.m_DT.is_infinite(add.at(i).v)
                || Vector_2(v_n.v->point(), q.v->point()).squared_length() <
                   Vector_2(add.at(i).v->point(), q.v->point()).squared_length())
                add.at(i) = v_n;   // if the saved vertex is infinite or longer than the current one, update
        }
        // p is guaranteed to be already added or will be added immediately after this step
        // so don't set isInZone to true until after it's passed
        if (v_n.key == p.key) isInZone = true;
    });

    for( SplitVertex v : add )
        if( !SG.m_DT.is_infinite( v.v ) ) {
            addPolygonSpannerEdge( E_P, q, v);
        }
}

inline void addPolygonEdges(SplitVertexEdgeMap& E_P,
                            SplitVertex& p,
                            SplitVertex& q,
                            SplitVertex& r ) {
    /*
     * Only add edges from P if they have not already been added.
     * Failure to check for this will result in breaking the
     * assumption that edges are only added to non-Complete
     * vertices.
     */
    if( !contains( E_P, p.key) || !contains( E_P.at(p.key), q.key ) )
        addPolygonSpannerEdge(E_P, p, q);
    if( !contains( E_P, q.key) || !contains( E_P.at(q.key), r.key ) )
        addPolygonSpannerEdge( E_P, q, r);
}

inline void addPolygonSpannerEdges(const DelaunayGraph& SG,
                                   SplitVertexSet& V,
                                   SplitVertexEdgeMap& E,
                                   SplitVertexEdgeMap& E_P,
                                   SplitVertex& p,
                                   SplitVertex& q,
                                   SplitVertex& r ) {
    //assert( !SG.m_DT.is_infinite(p.v) && !SG.m_DT.is_infinite(q.v) && !SG.m_DT.is_infinite(r.v) );
    if( p == r ) return;

    //cout<< " adding forward edges...\n";
    addForwardEdges(SG, V, E, E_P, p, q, r);
    //cout<< " adding cross edges...\n";
    addCrossEdges(SG, V, E, E_P, p, q, r);
}

void processVertex(const DelaunayGraph& SG,
                   SplitVertexSet& V,
                   SplitVertexEdgeMap& E,
                   SplitVertexEdgeMap& E_P,
                   SplitVertex& v_i) {
    SplitVertex s_1 = V.at(v_i.s_1),
                   s_m = get_s_m( SG, V, E, v_i ),
                   s_j,
                   s_k;

//        cout <<  " v_i:" << v_i
//             << " s_1:" << s_1
//             <<"\n";
    if( E_P.empty() ) {
//             cout<< " adding edges between s_1 and s_m...\n";
        addPolygonSpannerEdges(SG, V, E, E_P, s_1, v_i, s_m);
    } else {
        optional< index_t > s_j_key;
        index_t s_k_key = s_m.key;
        forEachNeighbor(SG, V, E_P, v_i, [&](SplitVertex &v_n) {
            if (!s_j_key) s_j_key = {v_n.key};
            s_k_key = v_n.key;
        });
        //assert(s_j_key);
        s_j =  V.at(*s_j_key );
        s_k =  V.at( s_k_key );
//            cout<< " s_j:" << s_j << " s_k:" << s_k <<"\n";
//            cout<< " adding edges between s_1 and s_j...\n";
        addPolygonSpannerEdges(SG, V, E, E_P, s_1, v_i, s_j);
//            cout<< " adding edges between s_k and s_m...\n";
        addPolygonSpannerEdges(SG, V, E, E_P, s_k, v_i, s_m);
    }

//        cout<< " s_m:" << s_m << "\n";
//        cout<< " adding edges from P...\n";
    addPolygonEdges( E_P, s_1, v_i, s_m);
}

} // namespace polygon_spanner

void PolygonSpanner( DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E ) {
    using namespace polygon_spanner;

    // Create a vertex status map
    VertexStatusMap status;

    queue< index_t > level; // BFS queue
    SplitVertexEdgeMap E_P;

    SplitVertex v_i = V.V.front();

    level.push( v_i.key );

    //SG.addToEventQueue( v_i, 0 );

    do { // loop through level queue
        v_i = V.at( level.front() );
//        cout<<"\nprocessing "<<v_i<<"\n";

        //if( !E_P.empty() ) assert( E_P.at( v_i.key ).size() <= 5 );

        processVertex(SG, V, E, E_P, v_i);

        // BFS housekeeping

        forEachNeighbor(SG, V, E, v_i,
            [&](const SplitVertex &v_n) {
                if (!contains(status, v_n.key)) { // If N[v_i] is NOT Known, queue it and add to Known
                    level.push(v_n.key);
                    status[v_n.key] = Known;
                }
            }
        );
        level.pop();
        status[v_i.key] = Complete;
    } while( !level.empty() ); // level is not empty

    VertexHandle v1;

    // Add all edges from E_P to SG
    for( auto& edge : E_P ) {
        // Each edge is a pair of size_t, unordered_set<size_t>
        v1 = V.at( edge.first ).v;
        for( auto v2 : edge.second ) {
            SG.addEdge(v1, V.at(v2).v);
        }
    }

    // Test degree assumption given after lemma 3.4
//    for( auto it=SG.m_E.begin(); it!=SG.m_E.end(); ++it ) {
//        assert( it->second.size() <= 27 );
//    }

    swap( E, E_P );

} // PolygonSpanner( SpanningGraph &P )

} // namespace bgs2005


template< typename RandomAccessIterator, typename OutputIterator >
void BGS2005( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result) {
    using namespace bgs2005;

    DelaunayGraph G( pointsBegin, pointsEnd ); // Step 1

    SpanningGraph( G ); // Step 2

    SplitVertexSet V;
    SplitVertexEdgeMap P;
    TransformPolygon( G, V, P ); // Step 3

    PolygonSpanner( G, V, P ); // Step 4

    // send resulting edge list to output iterator
    for( auto const& adj : G.m_E ) {
        VertexHandle v_1 = adj.first;
        for( auto const& v_2 : adj.second ) {
            *result = make_pair( v_1->info(), v_2->info() );
            ++result;
        }
    }
}

} // namespace spanner

#endif // SPANNERS_BGS2005_H
