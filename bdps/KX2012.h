
#ifndef SPANNERS_KX2012_H
#define SPANNERS_KX2012_H

#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

#include <boost/functional/hash.hpp> // size_t pair hash

#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include "printers/GraphPrinter.h"
#include "tools/DelaunayL2.h"
#include "tools/Metrics.h"
#include "tools/Utilities.h"


namespace spanner {

    using namespace std;

    namespace kx2012 {

        bool
        selectEdge(//const DelaunayL2 &T,
                   index_tPairMap &E,
                   const VertexHandle& i,
                   const VertexHandle& j) {
            //assert(T.is_edge(i, j));
            //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";
            auto edge_j_i = make_pair(j->info(), i->info());
            auto existing = E.begin();
            bool inserted = false;
            tie(existing, inserted) = E.try_emplace(make_pair(i->info(), j->info()), false);
            if (spanners::contains(E, edge_j_i)) { E[edge_j_i] = true; }

            return inserted;
        }

    } // namespace kx2012

    template<typename RandomAccessIterator, typename OutputIterator>
    void KX2012(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result,
                bool printLog = false) {
        using namespace kx2012;
        using spanners::contains;

        //if(printLog) cout<<"alpha:"<<alpha<<",";

        // Construct Delaunay triangulation

        vector<Point> P(pointsBegin, pointsEnd);
        vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayL2 T;

        //N is the number of vertices in the delaunay triangulation.
        index_t n = P.size();
        if (n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        vector<VertexHandle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        //cout<<"del:";
        {   //Timer tim;
            for (size_t entry : index) {
                auto vh = T.insert(P[entry], hint);
                hint = vh->face();
                vh->info() = entry;
                handles[entry] = vh;
            }
        }
        VertexHandle v_inf = T.infinite_vertex();
        index_tPairMap E; // list of potential edges, value must be true for insertion to result

        // Iterate through vertices in T

        for (auto m = T.finite_vertices_begin(); m != T.finite_vertices_end(); ++m) {
            //if( printLog ) cout<<"\n\nm:"<<m->info()<<" ";

            // Get neighbors of m
            VertexCirculator N = T.incident_vertices(m);

            if (T.is_infinite(N)) ++N; //Verify N is not starting at infinity
            VertexCirculator done(N); //Artificial end to circulator
            VertexHandle coneReferenceHandle;

            //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";



            //if(printLog) cout<<"done:"<<done->info()<<",";
            vector<VertexHandle> angleSet; //Store the 3 vertices in a vector for each getAngle set
            unordered_set<VertexHandle> wideVertices;
            Point p = m->point();
            index_t currentPointIndex = m->info();
            {
//            cout<<"a:";
//            Timer tim;
                while (angleSet.size() < 3) {
                    if (!T.is_infinite(N)) angleSet.emplace_back(N);
                    N++;
                }

                //N is set to
                number_t angleSum;

                if (T.is_infinite(N)) { N++; }
                coneReferenceHandle = N;
                do {
                    if (T.is_infinite(N)) { N++; }

                    angleSum = getAngle(angleSet[2]->point(), p, angleSet[0]->point());

                    if (angleSum > FOUR_PI_OVER_FIVE) {
                        wideVertices.insert(angleSet.begin(), angleSet.end());
                        if (printLog) {
                            cout << "Points: {" << angleSet[2]->info() << "," << m->info() << ","
                                 << angleSet[0]->info()
                                 << "} make getAngle: " << angleSum << endl;

                        }
                        coneReferenceHandle = angleSet[2];
                    }

                    angleSet[0] = angleSet[1];
                    angleSet[1] = angleSet[2];
                    angleSet[2] = N++;

                } while (angleSet[0] != done);


                for (const auto &v : wideVertices) {
                    selectEdge( E, m, v);
                }
            }
            Point coneReferencePoint(coneReferenceHandle->point());
            //number_t conalAngle;
            {
//            cout<<"b:";
//            Timer tim;
                if (T.is_infinite(N)) { N++; }
                done = N;
                unordered_map<cone_t, VertexHandle> closestVertexInCone;
                unordered_map<cone_t, number_t> closestPointDistanceInCone;

                do {
                    if (!T.is_infinite(N)) {
                        if (spanners::contains(wideVertices, N)) {
                            if (!closestVertexInCone.empty()) {
                                coneReferencePoint = N->point();
                                for (const auto &v : closestVertexInCone) {
                                    selectEdge( E, m, v.second);
                                }
                                closestVertexInCone.clear();
                                closestPointDistanceInCone.clear();
                            }
                        } else {
                            number_t conalAngle = getAngle(N->point(), p, coneReferencePoint);
                            cone_t currentCone = floor(conalAngle / PI_OVER_FIVE);
                            number_t currentDistance = getDistance(p, N->point());

                            if (!spanners::contains(closestPointDistanceInCone, currentCone) ||
                                (currentDistance < closestPointDistanceInCone[currentCone])) {
                                closestVertexInCone[currentCone] = N;
                                closestPointDistanceInCone[currentCone] = currentDistance;
                            }
                        }
                    }
                } while (++N != done);

                for (const auto &v : closestVertexInCone) {
                    selectEdge( E, m, v.second);
                }
            }

            {
//            cout<<"c:";
//            Timer tim;
                while (!spanners::contains(wideVertices, N) && ++N != done);
                done = N;
                pair<int, VertexHandle> previousPoint(-1, v_inf);
                do {
                    if (!T.is_infinite(N)) {

                        if (spanners::contains(wideVertices, N)) {
                            coneReferencePoint = N->point();
                            previousPoint.second = N;
                            previousPoint.first = -1;

                        } else {
                            number_t conalAngle = getAngle(N->point(), p, coneReferencePoint);
                            auto currentCone = cone_t(floor(conalAngle / PI_OVER_FIVE));
                            int conalDifference = int(currentCone) - int(previousPoint.first);

                            for (int conalEdgesToAdd = std::min(conalDifference - 1, 2);
                                 conalEdgesToAdd > 0; conalEdgesToAdd--) {
                                //Determine if one of the edges adjacent to the empty cone is selected already
                                bool containsPreviousPoint = spanners::contains(E,
                                                                                make_pair(currentPointIndex,
                                                                                        previousPoint.second->info()));
                                bool containsN = spanners::contains(E,
                                                                    make_pair(currentPointIndex, N->info()));
                                //assert(previousPoint.second != v_inf);
                                //assert(N->handle() != v_inf);
                                if (containsPreviousPoint && !containsN) {
                                    selectEdge(E, m, previousPoint.second);
                                } else if (!containsPreviousPoint && containsN) {
                                    selectEdge( E, m, N->handle());
                                } else //If neither adjacent edge is already selected, add the longest one
                                {
                                    VertexHandle edgeToAdd =
                                            getDistance(N->point(), p) > getDistance(previousPoint.second->point(), p)
                                            ? N->handle() : previousPoint.second;

                                    selectEdge( E, m, edgeToAdd);
                                }
                            }
                            previousPoint.first = int(currentCone);
                            previousPoint.second = N;
                        }
                    }
                } while (++N != done);
            }
        }

        // Done. Send edges from G_prime with value == true (selected by both endpoints) to output.

        // Edge list is only needed for printing. Remove for production.
//   vector< pair<size_t,size_t> > edgeList;
//   edgeList.reserve( E.size() );

        // Send resultant graph to output iterator
        for (auto e : E) {
            if (e.second) { // e.second holds the bool value of whether both vertices of an edge selected the edge
                // Edge list is only needed for printing. Remove for production.
//            edgeList.push_back(e.first);

                *result = e.first;
                ++result;
                //            *result = reverse_pair(e.first);
                //            ++result;
            }
        }


        //
        //
        // START PRINTER NONSENSE
        //
        //

//   if( printLog ) {
//       GraphPrinter printer(1);
//       GraphPrinter::OptionsList options;
//
//       options = {
//           { "color", printer.inactiveEdgeColor },
//           { "line width", to_string(printer.inactiveEdgeWidth) }
//       };
//       printer.drawEdges( T, options );
//
//       options = { // active edge options
//           { "color", printer.activeEdgeColor },
//           { "line width", to_string(printer.activeEdgeWidth) }
//       };
//       printer.drawEdges( edgeList.begin(), edgeList.end(),P, options );
//
//
//       options = {
//           { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
//           { "color", make_optional( printer.backgroundColor ) }, // text color
//           { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
//           { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
//       };
//       GraphPrinter::OptionsList borderOptions = {
//           { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
//           { "color", printer.activeEdgeColor }, // additional border color
//           { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
//       };
//       printer.drawVerticesWithInfo( T, options, borderOptions );
//
//       printer.print( "KX2012" );
//       cout<<"\n";
//   }

        //
        //
        // END PRINTER NONSENSE
        //
        //

    } // function KX2012

} // namespace spanner

#endif // SPANNERS_KX2012_H
