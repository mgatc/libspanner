#ifndef LIBSPANNER_BGHP2010_H
#define LIBSPANNER_BGHP2010_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <iostream>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//Boost library
//#include <boost/functional/hash.hpp> // size_t pair hash

//CGAL library
#include <CGAL/algorithm.h>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//Project library
#include "../constants.h"
#include "../delaunay/DelaunayTD.h"
#include "../bdps/types.h"
#include "../utilities.h"


namespace spanner {

    namespace bghp2010 {

        typedef DelaunayTD::VertexDescriptor VertexDescriptor;

        enum EdgeLabel { // per cone
            FIRST = 1,
            CLOSEST = 0,
            LAST = -1
        };

        inline number_t getCanonicalAngle(const index_t p, const index_t q, const index_t r, const bdps::input_t &P) {
            return CGAL::min(getAngle(p, q, r, P), getAngle(r, q, p, P));
        }

//Finds the bisector length of a given edge.
        inline number_t bisectorLength(const Edge &e, const std::vector<Point> &H) {
            cone_t cone = td::getCone(e.first, e.second, H);

//    //assert(cone<6);
//    //assert(e.first<H.size());

            number_t xCord = H[e.first].x() - td::orthBisectorSlopes[cone];
            number_t yCord = H[e.first].y() + 1;

            Point bisectorPoint(xCord, yCord);

            K::Line_2 bisectorLine(H[e.first], bisectorPoint);

            Point intersectionPoint = bisectorLine.projection(H[e.second]);

            number_t bisectorLen = getDistance(H[e.first], intersectionPoint);

            return bisectorLen;
        }

        inline cone_t iPlusOne(const cone_t i) {
            return ((i + 1) % 6);
        }

        inline cone_t iLessOne(const cone_t i) {
            return ((i - 1 + 6) % 6);
        }

        template<class KeyEdgesMap>
        void findKeyEdges(KeyEdgesMap &KeyEdges, DelaunayTD &D, const bdps::input_t &P) {
            typedef DelaunayTD::Edge_descriptor EdgeDescriptor;

            // find closest in each cone
            for (auto vit = D.finite_vertices_begin();
                 vit != D.finite_vertices_end(); ++vit) {
                auto w = *vit;

                for (size_t cone = 1; cone < 6; cone += 2) {
                    std::vector<EdgeDescriptor> fan;
                    D.fanOfCone(w, cone, fan);

                    // Find closest
                    VertexDescriptor closestKnown = w;
                    number_t closestKnownBisector = INF;

                    for (auto e : fan) {
                        auto v = D.source(e);
                        number_t length = bisectorLength(std::make_pair(w, v), P);
                        if (length < closestKnownBisector) {
                            closestKnown = v;
                            closestKnownBisector = length;
                        }
                    }
                    cone_t flattenedCone = cone / 2;
                    if (closestKnown != w) {
                        KeyEdges[CLOSEST][w][flattenedCone] = closestKnown;
                    }

                    // Find first and last
                    if (!fan.empty()) {
                        KeyEdges[FIRST][w][flattenedCone] = D.source(fan.back());
                        KeyEdges[LAST][w][flattenedCone] = D.source(fan.front());
                    }
                }
            }
        }

        template<class KeyEdgesContainer>
        bool iRelevant(const VertexDescriptor child,
                       const VertexDescriptor parent,
                       const VertexDescriptor grandparent,
                       cone_t i,
                       const bdps::input_t &P,
                       KeyEdgesContainer &KeyEdges) {
            if (child == SIZE_T_MAX || grandparent == SIZE_T_MAX)
                return false;

            cone_t childCone = td::getCone(child, grandparent, P);

            if (childCone != i)
                return false;

            return (child == KeyEdges[FIRST][parent][iPlusOne(i) / 2] &&
                    child != KeyEdges[CLOSEST][parent][iPlusOne(i) / 2])
                   || (child == KeyEdges[LAST][parent][iLessOne(i) / 2] &&
                       child != KeyEdges[CLOSEST][parent][iLessOne(i) / 2]);
        }

        template<class EdgeList, class KeyEdgesContainer>
        bool iDistant(const VertexDescriptor w,
                      const cone_t i,
                      const size_t charge,
                      const DelaunayTD &D,
                      const bdps::input_t &P,
                      const EdgeList &E,
                      KeyEdgesContainer &KeyEdges) {
            if (charge < 2)
                return false;

            VertexDescriptor u = D.parent(w, i),
                    first = KeyEdges[FIRST][w][iPlusOne(i) / 2],
                    last = KeyEdges[LAST][w][iLessOne(i) / 2];

            //cout<< "w:"<<w<<" u:"<<u<<" i:"<<i<< " first:"<<first<<" last:"<<last<<endl;

            //return false;
            return !contains(E, std::make_pair(w, u))
                   && iRelevant(first, w, u, iPlusOne(i), P, KeyEdges)
                   && iRelevant(last, w, u, iLessOne(i), P, KeyEdges);
        }

        template<class AdjacencyList, class ChargeList, class EdgeList>
        void addClosestInNegativeCones(const AdjacencyList &ClosestEdges,
                                       const bdps::input_t &P,
                                       ChargeList &Charges,
                                       EdgeList &E,
                                       bool printLog) {
            const index_t n = P.size();
            // Add closest in each negative cone for each vertex
            if (printLog) std::cout << "\nClosest\n";
            for (index_t u = 0; u < n; ++u) {
                for (auto v : ClosestEdges[u]) {
                    if (v != SIZE_T_MAX) {
                        cone_t cone = td::getCone(u, v, P);
                        auto pair1 = std::make_pair(u, cone),
                             pair2 = std::make_pair(v, (cone + 3) % 6);
                        Charges.try_emplace(pair1, 0);
                        Charges[pair1]++;
                        Charges.try_emplace(pair2, 0);
                        Charges[pair2]++;

                        E.emplace(u, v);
                        if (printLog) std::cout << v << " " << u << "\n";
                    }
                }
            }
        }

        template<class KeyEdgeList, class ChargeList, class EdgeList>
        void add_iRelevantNeighbors(KeyEdgeList &KeyEdges,
                                    const bdps::input_t &P,
                                    const DelaunayTD &D,
                                    ChargeList &Charges,
                                    EdgeList &E,
                                    bool printLog) {
            const index_t n = P.size();

            // Add first and last in each negative cone if it is (i+1)-relevant
            if (printLog) std::cout << "\nFirst and AlgorithmLast\n";
            for (index_t u = 0; u < n; ++u) {
                // get edges from positive cones
                for (auto it = D.positive_cone_edges_begin(u);
                     it != D.positive_cone_edges_end(u);
                     ++it) {
                    auto e = *it;
                    auto v = D.target(e);
                    cone_t posCone = td::getCone(u, v, P);
                    for (int j = -1; j <= 1; ++j) {
                        cone_t i = (posCone + 6 + j) % 6;
                        //size_t iLessOne = (i - 1 + 6) % 6;
                        auto w_label = static_cast<EdgeLabel>(j); // FIRST or LAST
                        auto w = KeyEdges[w_label][u][i / 2];

                        if (iRelevant(w, u, v, posCone, P, KeyEdges)) {

                            //size_t cone = getCone( u, v, P );
                            auto pair2 = std::make_pair(u, posCone);
                            Charges.try_emplace(pair2, 0);
                            Charges[pair2]++;

                            E.emplace(u, w);
                            if (printLog) std::cout << w << " " << u << "\n";
                        }
                    }
                }
            }
        }

        template<class KeyEdgeList, class ChargeList, class EdgeList>
        void handle_iDistantCharge2s(KeyEdgeList &KeyEdges,
                                     const bdps::input_t &P,
                                     const DelaunayTD &D,
                                     ChargeList &Charges,
                                     EdgeList &E) {
            //const size_t n = P.size();

            for (auto it = Charges.begin();
                 it != Charges.end(); ++it) {
                auto cone = *it;
                VertexDescriptor w = cone.first.first;
                cone_t i = cone.first.second;
                size_t &charge = cone.second;
                if (iDistant(w, i, charge, D, P, E, KeyEdges)) {
                    auto next = KeyEdges[FIRST][w][iPlusOne(i) / 2],
                            prev = KeyEdges[LAST][w][iLessOne(i) / 2];
                    E.emplace(next, prev);

                    // remove the edge to w that is farthest from w's grandparent in cone i
                    VertexDescriptor u = D.parent(w, i);
                    auto thetaNext = getCanonicalAngle(w, u, next, P),
                            thetaPrev = getCanonicalAngle(w, u, prev, P);
                    auto remove = std::make_pair(w, thetaNext > thetaPrev ? next : prev);
                    //assert(contains(E, remove));
                    E.erase(remove);

                    // update charges
                    auto nextCharge = std::make_pair(next, iPlusOne(i)),
                            prevCharge = std::make_pair(prev, iLessOne(i));
                    ++(Charges[nextCharge]);
                    ++(Charges[prevCharge]);
                    --charge;
                }
            }
        }

        template<class KeyEdgeList, class ChargeList, class EdgeList>
        void handleOtherCharge2s(KeyEdgeList &KeyEdges,
                                 const DelaunayTD &D,
                                 ChargeList &Charges,
                                 EdgeList &E,
                                 bool printLog) {
            //const size_t n = P.size();

            for (auto it = Charges.begin();
                 it != Charges.end(); ++it) {
                auto cone = *it;
                VertexDescriptor w = cone.first.first;
                cone_t i = cone.first.second;
                size_t &charge = cone.second;
                if (charge == 2
                    && Charges[std::make_pair(w, iLessOne(i))] == 1
                    && Charges[std::make_pair(w, iPlusOne(i))] == 1) {
                    auto w_parent = D.parent(w, i);
                    if (w_parent == SIZE_T_MAX)
                        continue;
                    auto next = KeyEdges[FIRST][w][iPlusOne(i) / 2],
                            prev = KeyEdges[LAST][w][iLessOne(i) / 2];

                    auto remove = std::make_pair(w, w == KeyEdges[LAST][w_parent][i] ? prev : next);

                    //assert(contains(E, remove));

                    if (printLog) std::cout << "Edge " << remove.first << "-" << remove.second << " ";
                    if (printLog && !contains(E, remove)) std::cout << "not ";
                    if (printLog) std::cout << "found\n";

                    E.erase(remove);

                    // update charge
                    --charge;
                }
            }
        }

    } // namespace bghp2010


// Main algorithm.
    void BGHP2010(const bdps::input_t& in, bdps::output_t& out,
                  bool printLog = false) {
        using namespace bghp2010;

        const index_t n = in.size();
        if (n > SIZE_T_MAX - 1 || n <= 1) return;

        // Step 1
        bdps::input_t P(in);
        DelaunayTD D(P.begin(), P.end());

        {
            //Timer tim;

            std::map<EdgeLabel, std::vector<std::vector<size_t>>> KeyEdges = {
                    {CLOSEST, std::vector<std::vector<index_t>>(n, std::vector<index_t>(3, SIZE_T_MAX))},
                    {FIRST,   std::vector<std::vector<index_t>>(n, std::vector<index_t>(3, SIZE_T_MAX))},
                    {LAST,    std::vector<std::vector<index_t>>(n, std::vector<index_t>(3, SIZE_T_MAX))}
            };

            findKeyEdges(KeyEdges, D, P);

            std::set<Edge> E;
            std::map<std::pair<VertexDescriptor, size_t>, size_t> Charges;

            // Step 2
            addClosestInNegativeCones(KeyEdges[CLOSEST], P, Charges, E, printLog);
            add_iRelevantNeighbors(KeyEdges, P, D, Charges, E, printLog);

            // Sanity check after step 2...
            //      each negative (odd) cone should have at most 1 edge,
            //      each positive (even) cone should have at most 2 edges

            if (printLog) std::cout << "\nCharges after step 2\n";
            for (auto charge: Charges) {
                cone_t cone = charge.first.second;
                //
                if (printLog)
                    std::cout << charge.first.first << "  "
                         << cone << "   "
                         << charge.second << "   "
                         << ((cone % 2 == 1 && charge.second <= 1)
                             || (cone % 2 == 0 && charge.second <= 2) ? "OK" : "FAIL")
                         << "\n";
                else
                    assert(charge.second <= 2 - (cone % 2));

            }

            // Step 3
            handle_iDistantCharge2s(KeyEdges, P, D, Charges, E);

            if (printLog) std::cout << "\nCharges after step 3\n";
            for (auto charge: Charges) {
                cone_t cone = charge.first.second;
                if (printLog) std::cout << charge.first.first << "  " << cone << "   " << charge.second << "   "
                                  << ((cone % 2 == 1 && charge.second <= 1) || (cone % 2 == 0 && charge.second <= 2)
                                      ? "OK" : "FAIL") << "\n";
            }
            //cout<<endl;

            // Step 4

            handleOtherCharge2s(KeyEdges, D, Charges, E, printLog);

            if (printLog) std::cout << "\nCharges after step 4\n";
            for (auto charge: Charges) {
                cone_t cone = charge.first.second;
                if (printLog) std::cout << charge.first.first << "  " << cone << "   " << charge.second << "   "
                                  << ((cone % 2 == 1 && charge.second <= 1) || (cone % 2 == 0 && charge.second <= 2)
                                      ? "OK" : "FAIL") << "\n";
            }

            // Send resultant graph to output iterator
            std::copy(E.begin(), E.end(), std::back_inserter(out));
//            for (auto e : E) {
//                *result = e;
//                ++result;
//            }





            // START PRINTER NONSENSE
//            if (printLog) {
//                vector<pair<Point, Point>> edgeList;
//
//                TikzPrinter printer("BGHP2010");
//                printer.autoscale(P.begin(), P.end(), 20);
//
//                TikzPrinter::OptionsList options;
//
//                options = {
//                        {"color",      printer.inactiveEdgeColor},
//                        {"line width", to_string(printer.inactiveEdgeWidth)}
//                };
//                printer.drawEdgesOfHalfTheta(D, options);
//
//                options = { // active edge options
//                        {"color",      printer.activeEdgeColor},
//                        {"line width", to_string(printer.activeEdgeWidth)}
//                };
//
//                printer.drawEdges(E.begin(), E.end(), P, options);
//
//                options = {
//                        {"vertex",     (to_string(printer.vertexRadius))}, // vertex width
//                        {"color",      (printer.backgroundColor)}, // text color
//                        {"fill",       (printer.activeVertexColor)}, // vertex color
//                        {"line width", (to_string(0))} // vertex border (same color as text)
//                };
//                TikzPrinter::OptionsList borderOptions = {
//                        {"border",     (to_string(printer.vertexRadius))}, // choose shape of vertex
//                        {"color",      printer.activeEdgeColor}, // additional border color
//                        {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
//                };
//                printer.drawVerticesWithInfo(D.points_begin(), D.points_end(), options, borderOptions);
//
//                printer.display();
//                cout << "\n";
//            }
            // END PRINTER NONSENSE
        }
    } // function BGHP2010

} // namespace spanner

#endif // SPANNERS_BGHP2010_H
