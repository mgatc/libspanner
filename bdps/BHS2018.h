//Needs optimizing currently testing.
#ifndef SPANNERS_BHS2017_H
#define SPANNERS_BHS2017_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//CGAL library
#include <CGAL/algorithm.h>
#include <CGAL/Line_2.h>

//Project library
#include "constants.h"
#include "delaunay/DelaunayL2.h"
#include "../bdps/types.h"
#include "Utilities.h"


namespace spanner {

    namespace bhs2018 {

        typedef CGAL::Line_2<K> Line;
        typedef std::unordered_map<Edge, number_t, IndexPairHash, IndexPairComparator> EdgeBisectorMap;
        typedef std::unordered_map<index_tPair, size_t, PointConeHash, PointConeComparator> PointConeMap;


        //Slopes of the cone boundry lines.
        //const vector<number_t> bisectorSlopes{INF, tan30, -1 * tan30, INF, tan30, -1 * tan30};
        const number_t orthBisectorSlopes[] = {0, -1 * COT30, COT30, 0, -1 * COT30, COT30};

        //Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an getAngle of (PI/3).
        inline cone_t getSingleCone(const index_t p, const index_t q, const vector<VertexHandle> &H) {
            const number_t alpha = PI / 3;
            const Point refPoint(H.at(p)->point().x() - TAN30, H[p]->point().y() + 1);
            //Point refPoint(H[p]->point().x(), H[p] ->point().y() + 1);

            number_t theta = getAngle(refPoint, H[p]->point(), H.at(q)->point());

            auto cone = cone_t(floor(theta / alpha));

            return cone;
        }

        //Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
        inline cone_t getCone(const index_t p, const index_t q, const vector<VertexHandle> &H) {
            return H[p] < H[q] ? getSingleCone(p, q, H) : (getSingleCone(q, p, H) + 3) % 6;
//            if (H[p] < H[q]) {
//                return getSingleCone(p, q, H);
//            } else {
//                return (getSingleCone(q, p, H) + 3) % 6;
//            }
        }

        //Finds the bisector length of a given edge.
        inline number_t bisectorLength(const Edge &e, const vector<VertexHandle> &H) {

            cone_t cone = getCone(e.first, e.second, H);
            //assert(cone < 6);
            //assert(e.first < H.size());

            number_t xCord = H[e.first]->point().x() - orthBisectorSlopes[cone];
            number_t yCord = H[e.first]->point().y() + 1;

            Point bisectorPoint(xCord, yCord);

            Line bisectorLine(H[e.first]->point(), bisectorPoint);

            Point intersectionPoint = bisectorLine.projection(H[e.second]->point());

            number_t bisectorLen = getDistance(H[e.first]->point(), intersectionPoint);

            return bisectorLen;
        }

        /*
          Step 3: Add incident edges consists of 2 sub steps.
          (3.1) Starts with the empty set E_A.
          (3.2) For each edge in L (sorted in non-decreasing order) Let i be the cone of p containing q if E_A has no edges with endpoint p in the
                neighborhood of p in cone i and E_A has no edges with endpoint q in the neighborhood of q in cone i+3 then add edge (p,q) to E_A.
        */
        inline void addIncident(std::vector<std::pair<size_t, size_t>> &E_A,
                                PointConeMap &AL_E_A,
                                const std::vector<VertexHandle> &h,
                                const std::vector<std::pair<index_tPair, number_t>> &l) {
            //Loops through the entire set L.
            for (auto e : l) {

                //Separates the edge (p,q) into the vertices p and q.
                const index_t &p = e.first.first;
                const index_t &q = e.first.second;

                //Computes the cone of p containing q.
                cone_t p_cone = getCone(p, q, h),
                        q_cone = getCone(q, p, h);

                /*Evaluates the emptiness of a cone, given a vertex in the set E_A.
                  If a cone is empty then the set E_A does not contain an edge with the
                  given endpoint in the cone calculated above, and the status will be
                  set to true.*/
                bool p_coneEmpty = AL_E_A.find(make_pair(p, p_cone)) == AL_E_A.end();
                bool q_coneEmpty = AL_E_A.find(make_pair(q, q_cone)) == AL_E_A.end();

                /*Checks that both cone neighborhood are empty, if these are both empty
                then the condition for step 3 is met and (p,q) is added to E_A. (3.2)*/
                if (p_coneEmpty && q_coneEmpty) {
                    E_A.push_back(e.first);

                    //Adds (p,q) to an adjacency list for future calculation.
                    AL_E_A.emplace(make_pair(p, p_cone), q);
                    AL_E_A.emplace(make_pair(q, q_cone), p);
                }
            }
        }

        inline void canonicalNeighborhood(std::vector<index_t> &canNeighbors,
                                          const index_t p,
                                          const index_t r,
                                          const cone_t cone,
                                          const DelaunayL2 &DT,
                                          const std::vector<VertexHandle> &H,
                                          const EdgeBisectorMap &B) {

            Edge e = make_pair(p, r);

            /*Vertex circulator oriented to r to find fist and last end vertex. Once r is the circulator is oriented to the first vertex in the cone,
              that is in the canonical neighborhood. Once found all neighbors are added in clockwise order. For a vertex to be in the canonical
              neighborhood it must have a bisector length greater than or equal to that of (p,r)
            */
            auto N_p = DT.incident_vertices(H[p]);

            while (++N_p != H[r]);

            while (!DT.is_infinite(++N_p) && getCone(p, N_p->info(), H) == cone &&
                   (B.at(make_pair(p, N_p->info())) > B.at(e)
                    || abs(B.at(make_pair(p, N_p->info())) - B.at(e)) < EPSILON));


            while (!DT.is_infinite(--N_p) && getCone(p, N_p->info(), H) == cone &&
                   (B.at(make_pair(p, N_p->info())) > B.at(e)
                    || abs(B.at(make_pair(p, N_p->info())) - B.at(e)) < EPSILON)) {
                canNeighbors.push_back(N_p->info());
            }
        }

        /*
          Step 4: Add canonical edges consists of 4 sub steps.
          (4.1) Let r be an element of the canonical neighborhood of p in the cone 0 of p.
          (4.2) Add inner edges if total neigborhood edges is 3 or more.
          (4.3) If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r.
          (4.4) Consider the first and last edge in the canoncial neighborhood. 3 criteria to add.
            (4.4 a) If the edges are in cone 1 or 5 with respect to a and z add.
            (4.4 b) If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add.
            (4.4 c) Checks if end edges have a end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
                    to a and z if found the edge (b,c) or (w,y) is added.
        */
        inline void addCanonical(std::vector<Edge> &E_CAN,
                                 const index_t p,
                                 const index_t r,
                                 const DelaunayL2 &DT,
                                 const std::vector<VertexHandle> &H,
                                 const EdgeBisectorMap &B,
                                 PointConeMap &AL_e_a) {

            //Creates an edge (p,r)
            //Edge e = make_pair(p, r);

            //Computes the cone of p containing r.
            cone_t p_cone = getCone(p, r, H);

            //Set of the canonical neighborhood of p in the cone of p containing r. (This cone will be considered as cone 0)
            std::vector<index_t> canNeighbors;

            canonicalNeighborhood(canNeighbors, p, r, p_cone, DT, H, B);
            //assert(!canNeighbors.empty());

            //Number of edges in the neighborhood.
            size_t canEdges = canNeighbors.size() - 1;

            //Must be at least 1 edge.
            if (canEdges > 1) {
                //Add inner edges if total neighborhood edges is 3 or more. (4.2)
                for (index_t i = 1; i < canEdges - 1; i++) {
                    E_CAN.emplace_back(canNeighbors.at(i), canNeighbors.at(i + 1));
                    //assert(DT.is_edge(H.at(canNeighbors.at(i)), H.at(canNeighbors.at(i + 1))));
                }

                //End edges in the canonical neighborhood.
                const std::vector<Edge> canExtrema{
                        make_pair(canNeighbors.at(1), canNeighbors.front()),
                        make_pair(canNeighbors.at(canEdges - 1), canNeighbors.back())
                };

                //If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r. (4.3
                for (auto edge : canExtrema) {
                    if (edge.second == r && canEdges > 1) {
                        E_CAN.push_back(edge);
                        //assert(DT.is_edge(H.at(edge.first), H.at(edge.second)));
                    }
                }

                //AlgorithmFirst and last edges in the canonical neighborhood are considered and added by 3 criteria. (4.4)
                std::vector<cone_t> cone(6);
                for (cone_t i = 0; i < 6; ++i) {
                    cone[i] = (p_cone + i) % 6;
                }

                //If the edges are in cone 1 or 5 with respect to a and z add. (4.4 a)
                for (size_t i = 0; i < canExtrema.size(); ++i) {
                    const auto edge = canExtrema[i];
                    const size_t z_cone = 1 + int(i == 1) * 4;
                    if (getCone(edge.second, edge.first, H) == cone[z_cone]) {
                        E_CAN.push_back(edge);
                        //if (printLog) cout << edge.first << "-" << edge.second << "[" << i << "]\\" << cone[z_cone] << "/,";

                        //assert(DT.is_edge(H.at(edge.second), H.at(edge.first)));
                    }
                }

                const std::vector<PointConeMap::iterator> endpointZ{
                        AL_e_a.find(std::make_pair(canExtrema[0].second, cone[2])),
                        AL_e_a.find(std::make_pair(canExtrema[1].second, cone[4]))
                };
                const auto blank = AL_e_a.end(); //Iterator to end of map to check if an edge exists.

                //If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add. (4.4 b)
                for (size_t i = 0; i < canExtrema.size(); ++i) {
                    const auto edge = canExtrema[i];
                    const cone_t z_cone = 2 + int(i == 1) * 2;
                    if (endpointZ[i] == blank && getCone(edge.second, edge.first, H) == cone[z_cone]) {
                        E_CAN.push_back(edge);

                        //assert(DT.is_edge(H.at(edge.first), H.at(edge.second)));
                    }
                }

                /*Checks if end edges have an end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
                  to a and z if found the edge (b,c) or (w,y) is added. (4.4 c)*/
                // (a,b)
                for (size_t i = 0; i < canExtrema.size(); ++i) {
                    const auto edge = canExtrema[i];
                    const cone_t z_cone = 2 + int(i == 1) * 2;
                    if (getCone(edge.second, edge.first, H) == cone[z_cone]
                        && endpointZ[i] != blank
                        && endpointZ[i]->second != edge.first) {
                        std::vector<size_t> zCanNeighbors;
                        canonicalNeighborhood(
                                zCanNeighbors, edge.second, endpointZ[i]->second,
                                cone[z_cone], DT, H, B
                        );
                        auto y = find(zCanNeighbors.begin(), zCanNeighbors.end(), edge.first);
                        auto w = y;
                        w += int(y == zCanNeighbors.begin());
                        w -= int(y == zCanNeighbors.end() - 1);

                        //assert(y != zCanNeighbors.end());
                        //assert(y != w);
                        //assert(DT.is_edge(H.at(*w), H.at(*y)));
                        E_CAN.emplace_back(*w, *y);
                    }
                }
            }
        }
    } // namespace BHS2018

//Main algorithm.
    void BHS2018(const bdps::input_t &in, bdps::output_t &out) {

        using namespace bhs2018;

        //Angle of the cones. Results in 6 cones for a given vertex.
        //const number_t alpha = PI / 3;

        bdps::input_t P(in);
        std::vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayL2 DT;

        //N is the number of vertices in the delaunay triangulation.
        size_t n = P.size();
        if (n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        std::vector<VertexHandle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        for (index_t entry : index) {
            auto vh = DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
            handles[entry] = vh;
        }

        //Put edges in a vector.
        std::vector<std::pair<Edge, number_t>> L;

        //Creates a map of edges as keys to its respective bisector length as the value. (Edges are not directional 1-2 is equivilent to 2-1)
        EdgeBisectorMap B(L.size());


        {
            //Timer t;

            for (auto e = DT.finite_edges_begin(); e != DT.finite_edges_end(); ++e) {
                auto edge = make_pair(
                        e->first->vertex((e->second + 1) % 3)->info(),
                        e->first->vertex((e->second + 2) % 3)->info()
                );
                auto length = bisectorLength(edge, handles);

                L.emplace_back(edge, length);
                B.emplace(edge, length);
            }
            //Timer t;
            //Step 2: Edges in the set L are sorted by their bisector length in non-decreasing order.
            std::sort(L.begin(), L.end(), [&](const auto &lhs, const auto &rhs) {
                return lhs.second < rhs.second;
            });
        }




        //Creates a set which will contain all edges returned by addIncident.
        std::vector<Edge> E_A;

        /*Creates an adjacency list where the inner lists are of size 6 representing the cones. The value stored in a particular inner index
          is the vertex that creates an edge with the outer vertex in the given cone. (i.e. If AL_E_A[10][4] = 5 in cone 4 of vertex 10 there
          exists an edge (10,5).)*/
        PointConeMap AL_E_A;

        //Step 3
        {
            //Timer t;
            addIncident(E_A, AL_E_A, handles, L);
        }

        //Add canonical E_CAN
        std::vector<Edge> E_CAN;

        //Step 4
        {

            //Timer t;
            for (auto e : E_A) {
                addCanonical(E_CAN, e.first, e.second, DT, handles, B, AL_E_A);
                addCanonical(E_CAN, e.second, e.first, DT, handles, B, AL_E_A);
            }
        }


        {
            //Timer t;
            //Union of sets E_A and E_CAN for final edge set removes duplicates.
            E_A.insert(E_A.end(), E_CAN.begin(), E_CAN.end());
            std::sort(E_A.begin(), E_A.end(), [](const auto &l, const auto &r) {
                return (CGAL::min(l.first, l.second) < min(r.first, r.second)
                        || (min(l.first, l.second) == min(r.first, r.second) &&
                            max(l.first, l.second) < max(r.first, r.second)));
            });
            E_A.erase(unique(E_A.begin(), E_A.end(), [](const auto &l, const auto &r) {
                return (l.first == r.first && l.second == r.second)
                       || (l.first == r.second && l.second == r.first);
            }), E_A.end());
        }

        // Edge list is only needed for printing. Remove for production.
        //vector<pair<Point,Point>> edgeList;

        // Send resultant graph to output iterator
        std::copy(E_A.begin(), E_A.end(), std::back_inserter(out));
//        for (auto e : E_A) {
//            // Edge list is only needed for printing. Remove for production.
//            //edgeList.emplace_back(handles.at(e.first)->point(), handles.at(e.second)->point());
//
//            *result = e;
//            ++result;
////        *result = make_pair(handles.at(e.second)->point(), handles.at(e.first)->point());
////        ++result;
//        }

        // START PRINTER NONSENSE
//    if(printLog) {
//        GraphPrinter printer(0.01);
//        GraphPrinter::OptionsList options;
//
//        options = {
//            {"color", printer.inactiveEdgeColor},
//            {"line width", to_string(printer.inactiveEdgeWidth)}
//        };
//        printer.drawEdges(DT, options);
//
//        options = { // active edge options
//            {"color", printer.activeEdgeColor},
//            {"line width", to_string(printer.activeEdgeWidth)}
//        };
//        printer.drawEdges(E_A.begin(), E_A.end(), P, options);
//
//
//        options = {
//            {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
//            {"color", make_optional(printer.backgroundColor)}, // text color
//            {"fill", make_optional(printer.activeVertexColor)}, // vertex color
//            {"line width", make_optional(to_string(0))} // vertex border (same color as text)
//        };
//        GraphPrinter::OptionsList borderOptions = {
//            {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
//            {"color", printer.activeEdgeColor}, // additional border color
//            {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
//        };
//        printer.drawVerticesWithInfo(DT, options, borderOptions);
//
//        printer.print("BHS2018");
//        cout << "\n";
//    }
        // END PRINTER NONSENSE

    } // function BHS2018

} // namespace spanner

#endif // SPANNERS_BHS2017_H
