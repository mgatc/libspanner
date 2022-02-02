#ifndef SPANNERS_BCC2012_H
#define SPANNERS_BCC2012_H

#include <bitset>
#include <cmath>         // ceil, floor, isinf
#include <functional>
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <utility>
#include <vector>        // handles

#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>


#include "constants.h"
#include "delaunay/DelaunayL2.h"
#include "../bdps/types.h"
#include "Utilities.h"



namespace spanner {

    namespace bcc2012 {

        const index_t TARGET = 7;

//        GraphPrinter tikz("../scratch/","scratch-graph");

        //set<Point> POINT_COLLECTOR;


//        GraphPrinter::OptionsList canonicalEdgeOptions = { // active edge options
//                {"color",      tikz.activeVertexColor},
//                {"line width", to_string(tikz.activeEdgeWidth)}
//        };

        struct WedgeParameters {
            index_t p;
            index_t q;
            cone_t cone;
        };

        std::pair<cone_t,cone_t> getCone(const std::vector<VertexHandle> &handles,
                                    const std::vector<index_t> &closest,
                                    const size_t p,
                                    const size_t q,
                                    const cone_t numCones) {
            const auto alpha = static_cast<number_t>(2 * PI / static_cast<number_t>(numCones));
            std::vector<cone_t> cone(3);
            for(int i=-1; i<2; ++i ) {
                cone[i+1] = static_cast<cone_t>((2*PI
                                                 + i*EPSILON
                                                 + getAngle(handles[closest[p]]->point(),
                                                           handles[p]->point(),
                                                           handles[q]->point()))
                                               / alpha) % numCones;
            }
            return make_pair( cone[1], cone[0]!=cone[1] ? cone[0] : cone[2] );
        }

//        inline cone_t getPreviousCone(const cone_t cone, const cone_t numCones) {
//            return (cone + numCones - 1) % numCones;
//        }

        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline bool vertexAgreesOnEdge(const std::vector<VertexHandle> &handles,
                                       std::vector<index_t> &closest,
                                       const std::vector<std::bitset<NUM_CONES>> &filled,
                                       const index_t p,
                                       const index_t q,
                                       std::pair<cone_t,cone_t> &cone) {
            if (closest.at(p) == SIZE_T_MAX) { // AlgorithmFirst, make sure the closest vertex is set
                closest[p] = q;
            }
            cone = getCone(handles, closest, p, q, NUM_CONES);

            bool pGivenConeFilled = filled.at(p)[cone.first],
                    pPrevConeFilled = filled.at(p)[cone.second];

            return !pGivenConeFilled && !pPrevConeFilled; ////////
        }

        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline void updateVertexConeStatus(std::vector<std::bitset<NUM_CONES>> &filled,
                                           std::vector<WedgeParameters> &wedge,
                                           const index_t p,
                                           const index_t q,
                                           const std::pair<cone_t,cone_t>& cone) {
            if (cone.second!=cone.first) {
                wedge.push_back({p, q, cone.second});
                filled.at(p)[cone.second] = true;
            }

            wedge.push_back({p, q, cone.first});
            filled.at(p)[cone.first] = true;
        }

        std::tuple<index_t,index_t,index_t> getNeighborsInCone(const index_t p, const index_t q, const cone_t cone, const cone_t numCones,
                                                          const DelaunayL2& DT, const std::vector<VertexHandle>& handles, const std::vector<index_t>& closest,
                                                               std::vector<index_t>& Q) {
            const auto ALPHA = static_cast<number_t>(2 * PI / static_cast<number_t>(numCones));

            // Line 2: Build Q
            auto N_p = DT.incident_vertices(handles[p]);

            while (++N_p != handles[q]); // point to q aka q_i

            while (!DT.is_infinite(++N_p) // move CCW until we leave the cone
                   && ((getCone(handles, closest, p, N_p->info(), numCones).first == cone
                        && getCone(handles, closest, p, N_p->info(), numCones).second == cone)
                       || N_p->info() == q));


            while (!DT.is_infinite(--N_p) // move CW until we leave the cone, adding each to Q
                   && ((getCone(handles, closest, p, N_p->info(), numCones).first == cone
                        && getCone(handles, closest, p, N_p->info(), numCones).second == cone)
                       || N_p->info() == q)) { // handles the case where q is on a boundary
                Q.push_back(N_p->info());
            }

            index_t j=0,
                    i=std::distance(Q.begin(),find(Q.begin(),Q.end(),q)),
                    k=Q.size()-1;

            return std::make_tuple(j,i,k);
        }
        void addEasyEdges(const index_t j, const index_t i, const index_t k,
                          const std::vector<index_t>& Q, std::vector<Edge>& add,
                          const std::unordered_set<index_t>& Forbidden = {}) {

            // Line 4: Add select edges
            if (i >= 2)
                for (size_t n = j+1; n <= i-2; ++n)
                    if (!contains(Forbidden, Q[n]) && !contains(Forbidden, Q[n+1]))
                        add.emplace_back(Q[n], Q[n+1]);

            if (k >= 2)
                for (size_t n = i+1; n <= k-2; ++n)
                    if (!contains(Forbidden, Q[n]) && !contains(Forbidden, Q[n+1]))
                        add.emplace_back(Q[n], Q[n+1]);
        }
/*
 *  Performs algorithm "wedge" with the exception of line 1. The
 *  for loop must be accomplished outside of the function.
 *  Degree 6 and 7 wedge algorithms are implemented as templates.
 *  This is the primary template, but will never be used. Instead,
 *  the specialized templates are used, defined below. clang-tidy
 *  complains about ununsed parameters here, it is safe to ignore.
 */
        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline void wedge(const DelaunayL2 &DT,
                          const std::vector<VertexHandle> &handles,
                          const std::vector<index_t> &closest,
                          const WedgeParameters &params,
                          std::vector<Edge> &addToE_star) {
            //assert(DEGREE == 6 || DEGREE == 7);
        }

        template<>
        inline void wedge<7>(const DelaunayL2 &DT,
                             const std::vector<VertexHandle> &handles,
                             const std::vector<index_t> &closest,
                             const WedgeParameters &params,
                             std::vector<Edge> &addToE_star) {
            const cone_t numCones = 8;

            const index_t &p = params.p;
            const index_t &q = params.q;
            const cone_t &cone = params.cone;
            index_t j,i,k;
            std::vector<index_t> Q;
            std::tie(j,i,k) = getNeighborsInCone(p,q,cone,numCones,DT,handles,closest,Q);

            std::vector<Edge> add;
            addEasyEdges(j,i,k,Q,add);
//            if (i>0) for (index_t m=j+1; m<i-1; ++m) {
//                addToE_star.emplace_back(Q[m],Q[m+1]);
//            }
//
//            if (k>0) for (index_t m=i+1; m<k-1 ;++m) {
//                addToE_star.emplace_back(Q[m],Q[m+1]);
//            }

            if (i!=j && i+1<k && getAngle(handles[Q[i+1]]->point(),handles[Q[i]]->point(),handles[p]->point()) > PI_OVER_TWO ) {
               add.emplace_back(Q[i],Q[i+1]);
            }
            if (i!=k && i>j+1 && getAngle(handles[p]->point(),handles[Q[i]]->point(),handles[Q[i-1]]->point()) > PI_OVER_TWO ) {
                add.emplace_back(Q[i],Q[i-1]);
            }

            for( auto e : add) {
                const index_t a = e.first,
                              b = e.second;
//                if( a==TARGET||b==TARGET ) {
//                    tikz.drawEdge(handles[a]->point(),handles[b]->point(), canonicalEdgeOptions);
//                    tikz.drawEdge(handles[p]->point(),handles[q]->point(), tikz.highlightEdgeOptions);
//                    tikz.drawCones(handles[p]->point(),handles[closest[p]]->point(), numCones, 1.0, tikz.coneOptions );
//                }
            }
            std::copy(add.begin(),add.end(),std::back_inserter(addToE_star));
        }

        template<>
        inline void wedge<6>(const DelaunayL2 &DT,
                             const std::vector<VertexHandle> &handles,
                             const std::vector<index_t> &closest,
                             const WedgeParameters &params,
                             std::vector<Edge> &addToE_star) {
            const cone_t numCones = 7;

            const index_t &p = params.p;
            const index_t &q = params.q;
            const cone_t &cone = params.cone;

            index_t j,i,k;
            std::vector<index_t> Q;
            std::tie(j,i,k) = getNeighborsInCone(p,q,cone,numCones,DT,handles,closest,Q);

//            if( find(Q.begin(),Q.end(),TARGET) != Q.end() ){
//                cout<<"Q contains the target!"<<"\n";
//            }

            // Line 3: Build Q'
            std::unordered_set<index_t> Q_prime; // select elements from Q (line 3)

            if( k>=3 ) { // if k is less than 3 then loop won't run.. if k is less than 1, k-1 overflows
                for( size_t m = j + 1; m < k - 1; ++m ) {
                    if( m == i )
                        continue;
                    auto theta = getAngle(handles[Q[m+1]]->point(),
                                          handles[Q[m]]->point(),
                                          handles[Q[m-1]]->point());
                    if( theta < SIX_PI_OVER_SEVEN) {
                        Q_prime.insert(Q[m]);
                    }
                }
            }

            std::vector<Edge> add;

            addEasyEdges(j,i,k,Q,add,Q_prime);

            for( auto e : add) {
                const index_t a = e.first,
                        b = e.second;
//                if( a==TARGET || b == TARGET ) {
//                    cout<<"\\texttt{Wedge}($"<<p<<","<<q<<"$) ($\\mathrm{cone}="<<cone<<", j="<<Q[j]<<", i="<<Q[i]<<", k="<<Q[k]<<"\\\\\n";
//                    for( size_t m=0; m<Q.size(); ++m ) {
//                        auto angle = getAngle(handles[Q[m+1]]->point(),
//                                              handles[Q[m]]->point(),
//                                              handles[Q[m-1]]->point());
//                        cout<<"$"<<Q[m]<<"$ & $"
//                            <<(0 < m && m < Q.size()-1 ? to_string(angle) : "\\mathrm{NA}") <<"$ & "
//                            <<(angle<SIX_PI_OVER_SEVEN&&m!=i&&m!=j&&m!=k ? "T":"F") << " & "
//                            <<(0<m &&( (Q[m]==b&&Q[m-1]==a) ||( Q[m]==a&&Q[m-1]==b)) ? "add" : "omit")<<"\\\\\n";
//                    }
//                    cout<<endl;
////                    tikz.drawEdge(handles[a]->point(),handles[b]->point(), canonicalEdgeOptions);
////                    tikz.drawEdge(handles[p]->point(),handles[q]->point(), tikz.highlightEdgeOptions);
//                    //tikz.drawCones(handles[p]->point(),handles[closest[p]]->point(), numCones, 1.0, tikz.coneOptions );
//
////                    vector<VertexHandle> pointsToAdd = { handles[p],handles[q],handles[a],handles[b] };
////                    transform(Q.begin(),Q.end(),inserter(POINT_COLLECTOR),[&](const auto& handle) {
////                        return handles[handle]->point();
////                    });
////                    POINT_COLLECTOR.insert(handles[p]->point());
////                    for(auto v : pointsToAdd ) {
////                        auto point = v->point();
////                        tikz.drawVertexWithLabel(point.x(),point.y(),to_string(v->info()), tikz.activeVertexOptions);
////                    }
//                }
            }

            std::copy(add.begin(),add.end(),back_inserter(addToE_star));
            // Line 5: if Q' is between j and i, reverse Q then apply the same procedure
            if(Q_prime.empty())
                return;

            add.clear();
            std::optional<size_t> iLessOneCcw = make_optional(p),
                             iLessOneCw =  i > 0 ? make_optional(Q[i-1]) : nullopt,
                             iPlusOneCcw = i < k ? make_optional(Q[i+1]) : nullopt,
                             iPlusOneCw = iLessOneCcw;

            std::function<bool(index_t)> inQPrime = [&](index_t v)-> bool {
                return contains(Q_prime,v);
            };
//            auto firstInQPrimeIterator = ;
            auto firstInQPrime = static_cast<size_t>(
                                 distance(Q.begin(),
                                          find_if(Q.begin(),Q.end(),inQPrime)) ),
                  lastInQPrime = static_cast<size_t>(
                                 distance(Q.begin(),
                                          find_if(Q.rbegin(),Q.rend(),inQPrime).base()) )
                                 - 1;

            assert (firstInQPrime > 0 && firstInQPrime < Q.size()-1);
            assert ( lastInQPrime > 0 &&  lastInQPrime < Q.size()-1);

            if( i > firstInQPrime ) {
                reverse(Q.begin(), Q.end());
                i = Q.size() - i - 1;

                std::swap(iLessOneCcw,iPlusOneCcw);
                std::swap(iLessOneCw, iPlusOneCw );
                //swap(j,k);
                std::swap(firstInQPrime,lastInQPrime);
            }

            // Line 6-7
            if(iLessOneCw && iLessOneCcw && i!=j && i-1 != j
               && getAngle(handles[*iLessOneCcw]->point(),
                           handles[q]->point(),
                           handles[*iLessOneCw]->point()) > FOUR_PI_OVER_SEVEN ) {
                add.emplace_back(Q[i],Q[i-1]);
            }

            // Line 8-9
            size_t f = firstInQPrime;
            size_t a = f+1;

            // Line 10-14
            if( iPlusOneCw && iPlusOneCcw && f == i+1 ) {
                number_t theta = getAngle(handles[*iPlusOneCcw]->point(),
                                          handles[q]->point(),
                                          handles[*iPlusOneCw]->point());
                if( a!=k
                    && ( theta < FOUR_PI_OVER_SEVEN || abs( theta - FOUR_PI_OVER_SEVEN) < EPSILON) ) {
                    add.emplace_back(Q[f],Q[a]);
                }
                if( f+1 != k && theta > FOUR_PI_OVER_SEVEN ) {
                    add.emplace_back(Q[i],Q[f+1]);
                }
            } else { // Lines 15-23
                size_t l = lastInQPrime,
                       b = l-1;

                if( l == k-1 ) {
                    add.emplace_back(Q[l],Q[b]);
                } else {
                    add.emplace_back(Q[b],Q[l+1]);
                    if( contains(Q_prime,Q[l-1]) ) {
                        add.emplace_back(Q[l],Q[l-1]);
                    }
                }
            }

//            for( auto e : add) {
//                const index_t a = e.first,
//                        b = e.second;
//                if( a==TARGET||b==TARGET ) {
//                    tikz.drawEdge(handles[a]->point(),handles[b]->point(), canonicalEdgeOptions);
//                    tikz.drawEdge(handles[p]->point(),handles[q]->point(), tikz.highlightEdgeOptions);
//                    tikz.drawCones(handles[p]->point(),handles[closest[p]]->point(), numCones, 1.0, tikz.coneOptions );
//
//                    vector<VertexHandle> pointsToAdd = { handles[p],handles[q],handles[a],handles[b] };
////                    for(auto v : pointsToAdd ) {
////                        auto point = v->point();
////                        tikz.drawVertexWithLabel(point.x(),point.y(),to_string(v->info()), tikz.activeVertexOptions);
////                    }
//                }
//            }
            std::copy(add.begin(),add.end(),back_inserter(addToE_star));
        }


    } // namespace bcc2012

    template<size_t DEGREE = 7, size_t NUM_CONES = DEGREE + 1>
    void BCC2012(const bdps::input_t& in, bdps::output_t& out) {
        using namespace bcc2012;



        //assert(DEGREE == 7 || DEGREE == 6);

        // Construct Delaunay triangulation
        bdps::input_t P(in);
        std::vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayL2 DT;

        //N is the number of vertices in the delaunay triangulation.
        index_t n = P.size();
        if (n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        std::vector<VertexHandle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        for (size_t entry : index) {
            auto vh = DT.insert(P[entry], hint);
            //cout<< entry<<" & "<< P[entry] << "\\\\\n";
            hint = vh->face();
            vh->info() = entry;
            handles[entry] = vh;
        }

        // Put edges in a vector, then sort on weight
        std::vector<Edge> L;

        for (auto e = DT.finite_edges_begin(); e != DT.finite_edges_end(); ++e) {
            L.emplace_back(e->first->vertex((e->second + 1) % 3)->info(),
                           e->first->vertex((e->second + 2) % 3)->info());
        }
        sort(L.begin(), L.end(), [&](const auto &lhs, const auto &rhs) {
            return getEdgeLength(lhs, P) < getEdgeLength(rhs, P);
        });

        //tikz.drawEdges(L.begin(),L.end(),P,tikz.triangulationEdgeOptions);

        std::vector<index_t> closest(n, SIZE_T_MAX); // id of closest vertex for orienting cones
        std::vector<std::bitset<NUM_CONES>> filled(n); // status of each cone
        std::vector<Edge> E; // output edge list
        std::vector<Edge> E_star; // edges added from "Wedge"

        for (auto pq : L) {
            index_t p = pq.first,
                   q = pq.second;

            // If either p or q's cone is filled, don't even bother
            if (filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES)
                continue;

            // Politely ask p if it wants an edge to q
            std::pair<cone_t,cone_t> cone_p(0,0);
            bool pAbides = vertexAgreesOnEdge<DEGREE>(handles, closest, filled, p, q,
                                                      cone_p);


            // Politely ask q if it wants an edge to p
            std::pair<cone_t,cone_t> cone_q(0,0);
            bool qAbides = vertexAgreesOnEdge<DEGREE>(handles, closest, filled, q, p,
                                                      cone_q);

//            cout<<p<<" & "
//                <<q<<" & "
//                <<getDistance(P[p],P[q])<<" & "
//
//                <<cone_p.first
//                <<(cone_p.first==cone_p.second?"":("/"+to_string(cone_p.second)))<<" & "
//                <<(filled[p][cone_p.first]?"F":"T")
//                <<(cone_p.first==cone_p.second ? "" : (filled[p][cone_p.first]?"/F":"/T"))<<" & "
//
//                <<cone_q.first
//                <<(cone_q.first==cone_q.second?"":("/"+to_string(cone_q.second)))<<" & "
//                <<(filled[p][cone_p.first]?"F":"T")
//                <<(cone_q.first==cone_q.second?"":(filled[q][cone_q.first]?"/F":"/T"))<<" & "
//
//                <<(pAbides&&qAbides? "add" : "omit")
//                <<"\\\\\n";
            // Only continue if p and q both consent to add the edge
            if (pAbides && qAbides) {
                E.emplace_back(p, q); // Place the edge
                std::vector<index_t> pqVec = {p, q};
                for(auto v : pqVec )
                    if(filled[v].count()==0) {
                        auto opposite = v == p ? q : p;
//                        tikz.drawCones(handles[v]->point(),
//                                       handles[opposite]->point(),
//                                       NUM_CONES, 1.0, tikz.coneOptions);
//                        tikz.drawEdge(handles[v]->point(),handles[opposite]->point(), tikz.activeEdgeOptions);
                    }
//                if( p==TARGET||q==TARGET ) {
//                    //tikz.drawEdge(handles[p]->point(), handles[q]->point(), tikz.activeEdgeOptions);
////                    if(filled[TARGET].count()==0)
////                        tikz.drawCones(handles[TARGET]->point(),
////                                       handles[p==TARGET?q:p]->point(),
////                                       NUM_CONES, 1.0, tikz.coneOptions);
////                    POINT_COLLECTOR.insert(handles[p]->point());
////                    POINT_COLLECTOR.insert(handles[q]->point());
//                    vector<VertexHandle> pointsToAdd = { handles[p],handles[q] };
////                    for(auto v : pointsToAdd ) {
////                        auto point = v->point();
////                        tikz.drawVertexWithLabel(point.x(),point.y(),to_string(v->info()), tikz.activeVertexOptions);
////                    }
//                }

                // Wedge on each cone of pqVec and qp
                // There will be at least one for each, but there could
                // be two cones for one or both pqVec and qp if the edgeat(0)
                // falls on the boundary of a cone and the cone is not already filled
                std::vector<WedgeParameters> W; // holds the parameters for each call to wedge

                // Bookkeeping for p
                updateVertexConeStatus<DEGREE>(filled, W, p, q, cone_p);

                // Bookkeeping for q
                updateVertexConeStatus<DEGREE>(filled, W, q, p, cone_q);

                std::vector<Edge> addToE_star;

                // Wedge on p, q
                for (auto params : W) {
                    // find q
                    auto q_z = DT.incident_vertices(handles.at(params.p));
                    while (++q_z != handles.at(params.q)); // point to q

                    wedge<DEGREE>(DT, handles, closest, params, addToE_star);
                }
                E_star.insert(E_star.end(), addToE_star.begin(), addToE_star.end());
            }
        }

        // Combine E and E_star, remove duplicates
        E.insert(E.end(), E_star.begin(), E_star.end());
        std::sort(E.begin(), E.end());
        E.erase(unique(E.begin(), E.end(), [](const auto &l, const auto &r) {
            return (l.first == r.first && l.second == r.second)
                   || (l.first == r.second && l.second == r.first);
        }), E.end());


//        GraphPrinter::OptionsList edgeOptions = { // active edge options
//                {"color",      tikz.activeEdgeColor},
//                {"line width", to_string(tikz.inactiveEdgeWidth/2)}
//        };
//        spanners::bcc2012::tikz.drawEdges(DT, tikz.triangulationEdgeOptions);
//
//         spanners::bcc2012::tikz.drawEdges(E.begin(),E.end(),P,edgeOptions);
//
//
//
////        vector<Point> target = {P[TARGET]};
//        spanners::bcc2012::tikz.drawVertices(P.begin(), P.end(), tikz.activeVertexOptions);
//
//


        // Send resultant graph to output iterator
        std::copy( E.begin(), E.end(), std::back_inserter(out) );

    } // function BCC2012



} // namespace spanner

#endif // SPANNERS_BCC2012_H
