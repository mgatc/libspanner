#ifndef LIBSPANNER_DEGREE3_H
#define LIBSPANNER_DEGREE3_H

#include <array>
#include <iostream>
#include <fstream>
#include <list>
#include <cassert>
#include <string>
#include <limits>
#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles
#include <optional>
#include <boost/functional/hash.hpp> // size_t pair hash : used in Yao_inf_4
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/squared_distance_2.h>
//#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

// includes for spatial sorting
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>


#include "constants.h"
#include "delaunay/DelaunayLinf.h"
#include "../bdps/types.h"
#include "Utilities.h"


namespace spanner {

    namespace degree3 {

        const number_t MAX_STRETCH_CONSTRAINT = 1000.0;
//
//        typedef size_t id_type;
//
//        template<typename T>
//        struct InfoConvert
//        {
//            typedef T Info;
//            typedef const Info& result_type;
//
//            inline const Info& operator()(const Info& info0, bool) const {
//                // just return the info of the supporting segment
//                return info0;
//            }
//            inline const Info& operator()(const Info& info0, const Info& , bool) const {
//                // just return the info of the supporting segment
//                return info0;
//            }
//        };
//        // functor that defines how to merge color info when a site (either
//        // point or segment) corresponds to point(s) on plane belonging to
//        // more than one input site
//        template<typename T>
//        struct InfoMerge
//        {
//            typedef T    Info;
//            typedef Info result_type;
//
//            inline Info operator()(const Info& info0, const Info& info1) const {
//                // just return the info of the supporting segment
//                return info0;
//            }
//        };
//
//        typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<K,CGAL::Field_with_sqrt_tag>  Gt;
//        typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
//                index_t,
//                InfoConvert<index_t>,
//                InfoMerge<index_t>> StorageTraits;
//        typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt, StorageTraits> LinfDelaunayGraph;
//        typedef LinfDelaunayGraph::Site_2 Site;
//        typedef Site::Point_2 Point;
        typedef DelaunayLinf::Vertex_circulator VertexCirculator;
        //typedef DelaunayLinf::Edge LinfEdge;
        typedef DelaunayLinf::Vertex_handle VertexHandle;
        typedef DelaunayLinf::Face_handle FaceHandle;

//        typedef CGAL::Spatial_sort_traits_adapter_2<K,
//                CGAL::Pointer_property_map<Point>::type> SearchTraits;


        enum AnchorType {
            None, Weak, Strong, StrongSelected, WeakSelected, StartOddChain
        };

        // project objects
        typedef std::vector<std::pair<VertexHandle, VertexHandle>> FanCones;
        typedef std::vector<std::pair<VertexHandle, index_t>> YaoCones;
        typedef std::vector<size_t> NumYaoEdges;
        typedef std::vector<std::pair<VertexHandle, AnchorType>> AnchorCones;
        typedef std::pair<VertexHandle, VertexHandle> SpannerEdge;
        typedef std::vector<std::vector<SpannerEdge>> SpannerCones;

        // inline functions required for bkpx2015

        cone_t getSingleCone(const VertexHandle &u, const VertexHandle &v) {
            const number_t alpha = PI / 2;
            const Point refPoint(u->site().point().x(), u->site().point().y() + 1);

            number_t theta = getAngle(refPoint, u->site().point(), v->site().point());
            auto cone = cone_t(floor(theta / alpha));

            return cone;
        }

        // get the cone of v wrt u (where u is at the center)
        cone_t getCone(const VertexHandle &u, const VertexHandle &v) {
            return u > v ? getSingleCone(u, v)
                         : (getSingleCone(v, u) + 2) % 4;
        }

        // add yao edges
        void addYaoEdges(std::vector<YaoCones> &yaoEdges,
                         std::vector<FanCones> &pointFans,
                         std::vector<NumYaoEdges> &yaoEdgeCount,
                         const std::vector<VertexHandle> &handles,
                         const DelaunayLinf &DT) {

            VertexCirculator circ = DT.incident_vertices(handles[0]); // default value
            std::vector<number_t> distances(4);

            cone_t cone = 0;
            index_t index = 0;

            for (const auto &point : handles) {

                circ = DT.incident_vertices(point);
                distances = {INF, INF, INF, INF};
                index = point->storage_site().info();

                YaoCones &edges = yaoEdges[index];
                FanCones &fans = pointFans[index];

                while (DT.is_infinite(circ))
                    ++circ;

                cone_t previousCone = getCone(point, circ);

                auto startPoint = circ++;

                while (!DT.is_infinite(circ) && getCone(point, circ) == previousCone && circ != startPoint)
                    ++circ;

                while (DT.is_infinite(circ))
                    ++circ;

                auto endpoint = circ;
                cone = getCone(point, circ);
                previousCone = 10; // invalid value, default as 10

                index_t point_id = point->storage_site().info();

                do {

                    if (!DT.is_infinite(circ)) {
                        cone = getCone(point, circ);
                        ++(edges[cone].second);

                        if (cone != previousCone)
                            fans[cone].first = circ;

                        fans[cone].second = circ;
                        previousCone = cone;

                        number_t proposedDistance = CGAL::l_infinity_distance(point->site().point(), circ->site().point());

                        bool tie = false;
                        if (abs(proposedDistance - distances[cone]) < EPSILON) {
                            tie = true;
                            auto v = edges[cone].first;
                            index_t v_id = v->storage_site().info();
                            index_t circ_id = circ->storage_site().info();

                            if (circ_id < v_id) {
                                distances[cone] = proposedDistance;
                                edges[cone].first = circ;
                            }
                        }

                        if (proposedDistance < distances[cone] && !tie) {
                            distances[cone] = proposedDistance;
                            edges[cone].first = circ;
                        }
                    }

                } while (++circ != endpoint); // finished determining the Yao edges + how many points are in fan of u's cone i
//            }
//
//            for (const auto &point : handles) {

                for (cone_t otherCone = 0; otherCone < 4; otherCone++) {

                    if (yaoEdges[index][otherCone].second == 0)
                        continue;

                    auto v = yaoEdges[index][otherCone].first;
                    index_t v_id = v->storage_site().info();

                    if (yaoEdges[v_id][(otherCone + 2) % 4].first != point) {
                        ++(yaoEdgeCount[index][otherCone]);
                        ++(yaoEdgeCount[v_id][(otherCone + 2) % 4]);
                    }
                }
            }
        }


        void determineAnchors(std::vector<AnchorCones> &anchorEdges,
                              std::vector<YaoCones> &yaoEdges,
                              std::vector<FanCones> &pointFans,
                              std::vector<NumYaoEdges> &yaoEdgeCount,
                              std::vector<VertexHandle> &handles,
                              DelaunayLinf &DT) {

            for (const auto &u : handles) {

                index_t u_id = u->storage_site().info();

                AnchorCones &anchors = anchorEdges[u_id];

                for (cone_t cone = 0; cone < 4; cone++) {

                    // no neighbors in cone
                    if (yaoEdges[u_id][cone].second == 0) {
                        continue;
                    }

                    auto v = yaoEdges[u_id][cone].first;

                    // only one neighbor in cone --> check if mutually single
                    if (yaoEdgeCount[u_id][cone] == 1) {
                        if (yaoEdgeCount[v->storage_site().info()][(cone + 2) % 4] == 1) {
                            anchorEdges[u_id][cone].first = v;
                            anchorEdges[u_id][cone].second = Weak;
                        }
                    }

                    if (yaoEdgeCount[u_id][cone] >= 2) {

                        auto v1 = pointFans[u_id][cone].first;
                        auto vk = pointFans[u_id][cone].second;
                        VertexCirculator current = DT.incident_vertices(u);
                        while (current != v1) { ++current; }
                        size_t numNeighbors = yaoEdges[u_id][cone].second;
                        size_t lcount = 1;

                        while (current != v) {
                            ++lcount;
                            ++current;
                        }

                        // establish vlower and vhigher and recalibrate circ accordingly
                        auto vlower = current;
                        --vlower;
                        auto vhigher = current;
                        ++vhigher;

                        bool ccwCanonical = (lcount >= 2 && !(DT.is_infinite(vlower)) &&
                                             (yaoEdges[vlower->storage_site().info()][getCone(vlower, v)].first == v &&
                                              yaoEdges[v->storage_site().info()][getCone(v, vlower)].first != vlower));

                        bool cwCanonical = (!ccwCanonical && lcount <= (numNeighbors - 1) &&
                                            !(DT.is_infinite(vhigher)) &&
                                            (yaoEdges[vhigher->storage_site().info()][getCone(vhigher, v)].first == v) &&
                                            yaoEdges[v->storage_site().info()][getCone(v, vhigher)].first != vhigher);

                        bool inCanonical = ccwCanonical || cwCanonical;
                        size_t position = lcount;
                        auto previous = current;
                        auto crown = previous;

                        int direction = 1 - 2 * ccwCanonical;

                        if (inCanonical) {
                            while (inCanonical) {
                                previous = current;

                                if (ccwCanonical)
                                    --current;
                                else
                                    ++current;

                                position += direction;

                                bool inFan = !((ccwCanonical && position < 1) ||
                                               (cwCanonical && position > numNeighbors));
                                bool unidirectional = (!(DT.is_infinite(current)) &&
                                                       yaoEdges[current->storage_site().info()][getCone(current,previous)].first == previous &&
                                                       yaoEdges[previous->storage_site().info()][getCone(previous,current)].first != current);

                                inCanonical = inFan && unidirectional;

                                bool yaoConnected = inCanonical && (yaoEdges[u_id][cone].first == current ||
                                                                    yaoEdges[current->storage_site().info()][(cone + 2) % 4].first == u);

                                if (yaoConnected)
                                    crown = current;

                            }

                            assert(yaoEdges[u_id][cone].first == crown || yaoEdges[crown->storage_site().info()][(cone + 2) % 4].first == u);

                            anchorEdges[u_id][cone].first = crown;
                            anchorEdges[u_id][cone].second = Weak;

                        } else {
                            anchorEdges[u_id][cone].first = v;
                            anchorEdges[u_id][cone].second = Weak;
                        }
                    } // when there are more than 2 edges!
                }
            }

            for (const auto &u : handles) {
                for (cone_t cone = 0; cone < 4; cone++) {
                    index_t u_id = u->storage_site().info();
                    VertexHandle v = anchorEdges[u_id][cone].first;

                    if (!DT.is_infinite(v) && (anchorEdges[v->storage_site().info()][(cone + 2) % 4].first == u
                                               || DT.is_infinite(anchorEdges[v->storage_site().info()][(cone + 2) % 4].first))) {
                        anchorEdges[u_id][cone].second = Strong;
                    }
                }

            } // strong anchors are identified


            // now it is time to select anchors
            for (const auto &u : handles) {
                index_t u_id = u->storage_site().info();

                for (cone_t cone = 0; cone < 4; cone++) {
                    if (DT.is_infinite(anchorEdges[u_id][cone].first))
                        continue;

                    if (anchorEdges[u_id][cone].second == Strong)
                        anchorEdges[u_id][cone].second = StrongSelected;

                    if (anchorEdges[u_id][cone].second == Weak) {
                        VertexCirculator circ = DT.incident_vertices(u);
                        auto v1 = pointFans[u_id][cone].first,
                                vk = pointFans[u_id][cone].second;
                        while (circ != v1) { ++circ; }

                        bool found = false;

                        do {
                            found = anchorEdges[circ->storage_site().info()][(cone + 2) % 4].first == u &&
                                    (anchorEdges[circ->storage_site().info()][(cone + 2) % 4].second == Weak ||
                                     anchorEdges[circ->storage_site().info()][(cone + 2) % 4].second == WeakSelected);
                        } while (!found && circ++ != vk);

                        if (found) { continue; }

                        std::vector<VertexHandle> visited;
                        auto previous = u;
                        auto current = anchorEdges[u_id][cone].first;
                        index_t previous_id = u_id;
                        cone_t currentCone = cone;
                        cone_t localCone = (cone + 2) % 4;
                        bool inChain = true;

                        do {

                            visited.push_back(previous);
                            previous = current;
                            previous_id = previous->storage_site().info();

                            currentCone = (visited.size() % 2) ? localCone : cone;

                            inChain = !(anchorEdges[previous_id][currentCone].second == Strong || anchorEdges[previous_id][currentCone].second == StrongSelected);

                            if (!inChain)
                                break;

                            current = anchorEdges[previous_id][currentCone].first;

                            assert(!DT.is_infinite(current));

                        } while (inChain);

                        assert(anchorEdges[previous_id][currentCone].second == Strong || anchorEdges[previous_id][currentCone].second == StrongSelected);

                        bool oddChain = (visited.size() % 2);

                        size_t position = 0;
                        current = visited.at(position);
                        index_t current_id = current->storage_site().info();

                        currentCone = cone;

                        if (oddChain) {
                            anchorEdges[current_id][cone].second = StartOddChain;
                            currentCone = localCone;
                            ++position;
                        }

                        for (size_t i = position; i < visited.size(); i += 2) {
                            current = visited.at(i);
                            current_id = current->storage_site().info();
                            anchorEdges[current_id][currentCone].second = WeakSelected;
                        }
                    }
                }
            }
        } // function Complete

        bool inEdgeList(const std::vector<SpannerEdge> &edgeList,
                        const VertexHandle u,
                        const VertexHandle v) {

            return any_of(edgeList.begin(), edgeList.end(), [&](const auto &e) {
                return (u == e.first || u == e.second) && (v == e.first || v == e.second);
            });
        }


        void degreeEightSpanner(std::vector<SpannerCones> &H8,
                                std::vector<AnchorCones> &anchorEdges,
                                std::vector<YaoCones> &yaoEdges,
                                std::vector<FanCones> &pointFans,
                                std::vector<NumYaoEdges> &yaoEdgeCount,
                                std::vector<VertexHandle> &handles,
                                DelaunayLinf &DT) {
            // put the edges conal vector
            for (const auto &w : handles) {
                SpannerCones edges(4);
                H8[w->storage_site().info()] = edges;
            }

            for (const auto &w : handles) {

                index_t w_id = w->storage_site().info();

                for (cone_t cone = 0; cone < 4; cone++) {

                    if (yaoEdges[w_id][cone].second == 0) { continue; }

                    if (anchorEdges[w_id][cone].second == StrongSelected ||
                        anchorEdges[w_id][cone].second == WeakSelected) {

                        auto v = anchorEdges[w_id][cone].first;
                        bool original = !inEdgeList(H8[w_id][cone], w, v);

                        if (original) {
                            SpannerEdge edge = std::make_pair(w, v);
                            H8[w_id][cone].push_back(edge);
                            H8[v->storage_site().info()][(cone + 2) % 4].push_back(edge);
                        }

                    }

                    if (yaoEdges[w_id][cone].second > 1) {

                        VertexCirculator current = DT.incident_vertices(w);
                        while (current != pointFans[w_id][cone].first) { ++current; }
                        auto previous = current;
                        ++current;

                        size_t position = 1;
                        size_t total = yaoEdges[w_id][cone].second;

                        while (position < total) {

                            cone_t cwCone = getCone(current, previous);
                            cone_t ccwCone = getCone(previous, current);

                            index_t current_id = current->storage_site().info();
                            index_t previous_id = previous->storage_site().info();

                            bool yaoConnected = (yaoEdges[w_id][cone].first == current ||
                                                 yaoEdges[current_id][(cone + 2) % 4].first == w)
                                                && (yaoEdges[w_id][cone].first == previous ||
                                                    yaoEdges[previous_id][(cone + 2) % 4].first == w);

                            bool uniCW = yaoEdges[current_id][cwCone].first == previous
                                         && yaoEdges[previous_id][ccwCone].first != current;

                            bool uniCCW = yaoEdges[current_id][cwCone].first != previous
                                          && yaoEdges[previous_id][ccwCone].first == current;

                            bool nonAnchor = anchorEdges[current_id][cwCone].first != previous
                                             && anchorEdges[previous_id][ccwCone].first != current;

                            bool addEdge = (uniCW || uniCCW) && yaoConnected;

                            if (yaoConnected && position == 1 && uniCW && nonAnchor) {

                                bool dual =
                                        yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[previous_id][(cone + 2) % 4] > 1;

                                bool startOdd = anchorEdges[previous_id][(cone + 2) % 4].first == w
                                                && anchorEdges[previous_id][(cone + 2) % 4].second == StartOddChain;

                                addEdge = !(dual && !startOdd);

                            }

                            if (yaoConnected && position == total - 1 && uniCCW && nonAnchor) {

                                bool dual =
                                        yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[current_id][(cone + 2) % 4] > 1;

                                bool startOdd = anchorEdges[current_id][(cone + 2) % 4].first == w
                                                && anchorEdges[current_id][(cone + 2) % 4].second == StartOddChain;

                                addEdge = !(dual && !startOdd);

                            }

                            if (addEdge) {

                                auto source = previous;
                                auto target = current;
                                size_t directedCone = ccwCone;

                                if (uniCW) {

                                    source = current;
                                    target = previous;
                                    directedCone = cwCone;

                                }

                                index_t source_id = source->storage_site().info();
                                index_t target_id = target->storage_site().info();

                                if (!inEdgeList(H8[source_id][directedCone], source, target)) {

                                    H8[source_id][directedCone].emplace_back(source, target);

                                    if (nonAnchor) {
                                        H8[target_id][(cone + 2) % 4].emplace_back(source, target);
                                    } else {
                                        H8[target_id][(directedCone + 2) % 4].emplace_back(source, target);
                                    }

                                }

                            }

                            ++position;
                            previous = current;
                            ++current;

                        }
                    }
                }
            }
        }


        void processSpanner(std::vector<SpannerCones> &H8,
                            const std::vector<AnchorCones> &anchorEdges,
                            const std::vector<YaoCones> &yaoEdges,
                            const std::vector<FanCones> &pointFans,
                            const std::vector<VertexHandle> &handles,
                            const DelaunayLinf &DT) {

            for (const auto &u : handles) {

                index_t u_id = u->storage_site().info();

                for (cone_t cone = 0; cone < 4; cone++) {

                    size_t charge = H8[u_id][cone].size();

                    if (charge == 1) {

                        SpannerEdge edge = H8[u_id][cone][0];

                        auto source = edge.first;
                        auto target = edge.second;

                        index_t source_id = source->storage_site().info();
                        index_t target_id = target->storage_site().info();

                        bool unidirectional = yaoEdges[source_id][cone].first == target &&
                                              yaoEdges[target_id][(cone + 2) % 4].first != source;
                        bool nonanchor = anchorEdges[source_id][cone].first != target &&
                                         anchorEdges[target_id][(cone + 2) % 4].first != source;

                        if (!(unidirectional && nonanchor)) { continue; }

//                    cout << endl << endl << "duplicate edge chain begins at " << source_id << endl;

                        for (cone_t i = 1; i <= 3; i += 2) {

                            cone_t localCone = (cone + i) % 4;

                            if (inEdgeList(H8[target_id][localCone], source, target) &&
                                H8[target_id][localCone].size() == 2) {

                                std::vector<VertexHandle> visited;
                                auto previous = source;
                                index_t previous_id = previous->storage_site().info();
                                auto current = target;
                                index_t current_id = current->storage_site().info();
                                cone_t currentCone = cone;
                                cone_t indicatorCone = localCone;

                                do {

                                    visited.push_back(previous);

                                    previous = current;
                                    previous_id = previous->storage_site().info();

                                    currentCone = (visited.size() % 2) ? localCone : cone;
                                    current = yaoEdges[previous_id][currentCone].first;

                                } while (inEdgeList(H8[previous_id][currentCone], previous, current) &&
                                         H8[previous_id][currentCone].size() == 2);

                                size_t total = visited.size() - 1;

                                // it's time to process edges

                                for (size_t j = 0; j < total; j += 2) {

                                    size_t l = (total - j);

                                    auto finish = visited.at(l);
                                    auto start = visited.at(l - 1);

                                    index_t finish_id = finish->storage_site().info();
                                    index_t start_id = start->storage_site().info();

                                    SpannerEdge newEdge = std::make_pair(start, finish);

                                    for (cone_t k = 0; k < 4; k++) {

                                        while (inEdgeList(H8[start_id][k], start, finish))
                                            H8[start_id][k].erase(
                                                    find(H8[start_id][k].begin(), H8[start_id][k].end(), newEdge));

                                        while (inEdgeList(H8[finish_id][k], start, finish))
                                            H8[finish_id][k].erase(
                                                    find(H8[finish_id][k].begin(), H8[finish_id][k].end(), newEdge));

                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (const auto &u : handles) {

                index_t u_id = u->storage_site().info();
                VertexCirculator middle = DT.incident_vertices(u);

                for (cone_t cone = 0; cone < 4; cone++) {

                    if (yaoEdges[u_id][cone].second >= 3) {

                        auto v1 = pointFans[u_id][cone].first,
                                vk = pointFans[u_id][cone].second;

                        while (middle != v1) { ++middle; }
                        ++middle;

                        // will continue so long as the charge is not 2 on middle in cone i+2 and middle is not vk

                        size_t total = yaoEdges[u_id][cone].second,
                                position = 1;

                        bool found = false;

                        while (!found && position < total - 1) {

                            if (H8[middle->storage_site().info()][(cone + 2) % 4].size() == 2 && middle != vk) {
                                auto previous = middle,
                                        next = middle;

                                --previous;
                                ++next;

                                index_t previous_id = previous->storage_site().info();
                                index_t middle_id = middle->storage_site().info();
                                index_t next_id = next->storage_site().info();

                                cone_t previousCone = getCone(previous, middle);
                                cone_t nextCone = getCone(next, middle);

                                bool edgePair = (yaoEdges[previous_id][previousCone].first == middle &&
                                                 yaoEdges[middle_id][(previousCone + 2) % 4].first != previous
                                                 && anchorEdges[previous_id][previousCone].first != middle &&
                                                 anchorEdges[middle_id][(previousCone + 2) % 4].first != previous)
                                                && (yaoEdges[next_id][nextCone].first == middle &&
                                                    yaoEdges[middle_id][(nextCone + 2) % 4].first != next
                                                    && anchorEdges[next_id][nextCone].first != middle &&
                                                    anchorEdges[middle_id][(nextCone + 2) % 4].first != next);

                                if (edgePair) {

                                    SpannerEdge previousEdge = std::make_pair(previous, middle);
                                    SpannerEdge nextEdge = std::make_pair(next, middle);

                                    //                         remove previousEdge
                                    while (inEdgeList(H8[previous_id][previousCone], previous, middle)) {
                                        H8[previous_id][previousCone].erase(find(H8[previous_id][previousCone].begin(),
                                                                                 H8[previous_id][previousCone].end(),
                                                                                 previousEdge));
                                    }


                                    while (inEdgeList(H8[middle_id][(cone + 2) % 4], previous, middle)) {
                                        H8[middle_id][(cone + 2) % 4].erase(find(H8[middle_id][(cone + 2) % 4].begin(),
                                                                                 H8[middle_id][(cone + 2) % 4].end(),
                                                                                 previousEdge));
                                    }

                                    // remove nextEdge
                                    while (inEdgeList(H8[next_id][nextCone], next, middle)) {
                                        H8[next_id][nextCone].erase(
                                                find(H8[next_id][nextCone].begin(), H8[next_id][nextCone].end(),
                                                     nextEdge));
                                    }

                                    while (inEdgeList(H8[middle_id][(cone + 2) % 4], next, middle)) {
                                        H8[middle_id][(cone + 2) % 4].erase(find(H8[middle_id][(cone + 2) % 4].begin(),
                                                                                 H8[middle_id][(cone + 2) % 4].end(),
                                                                                 nextEdge));
                                    }

                                    //                        cout << "adding shortcut: <" << previous_id << ", " << next_id << "> " << endl;

                                    SpannerEdge shortcut = std::make_pair(previous, next);
                                    H8[previous_id][previousCone].push_back(shortcut);
                                    H8[next_id][nextCone].push_back(shortcut);

                                }
                            }

                            ++position;
                            ++middle;

                        }
                    }
                }
            }
        }

//    bool BFS(vector<spannerCones> &H8, const Vertex_handle u, const Vertex_handle v) {
//
//        unordered_set<Vertex_handle> visited;
//        queue<Vertex_handle> next;
//
//        auto goal = v;
//        size_t waveTotal = 1;
//        size_t depth = 0;
//        next.push(u);
//
//        while (depth <= 50) { // depth limited BFS
//
//            size_t waveCount = 0;
//            size_t nextTotal = 0;
//
//            // process each wave
//            while (waveCount < waveTotal) {
//
//                Vertex_handle current = next.front();
//                next.pop();
//                size_t current_id = current->storage_site().info();
//
//                // if current has been visited continue
//                if (visited.count(current)) {
//                    waveCount++;
//                    continue;
//                }
//
//                visited.insert(current);
//
//                // add all unvisited vertices neighboring current in this wave to next
//                for (size_t cone = 0; cone < 4; cone++) {
//
//                    for (auto edge : H8[current_id][cone]) {
//
//                        auto neighbor = (edge.first == current) ? edge.second : edge.first;
//
//                        if (neighbor == goal)
//                            return true;
//
//                        if (visited.count(neighbor))
//                            continue;
//
//                        nextTotal++;
//                        next.push(neighbor);
//                    }
//                }
//
//                waveCount++;
//            }
//
//            waveTotal = nextTotal;
//            ++depth;
//
//        }
//
//        return false;
//
//    }

        template< typename Point >
        struct EuclideanDistanceToPoint {
            Point goal;
            explicit EuclideanDistanceToPoint(const Point& p) : goal(p) {}
            number_t operator()( Point p ) {
                return (number_t) sqrt(CGAL::squared_distance(p,goal));
            }
        };

        template<typename T>
        struct MinHeapCompare {
            bool operator()( const T &n1, const T &n2 ) const {
                return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
            }
        };

//        template<typename T>
//        struct MaxHeapCompare {
//            bool operator()( const T &n1, const T &n2 ) const {
//                return (n1.first < n2.first) || ((n1.first == n2.first) && (n1.second < n2.second));
//            }
//        };

        template< typename VertexContainer >
        std::optional<number_t> AStar( VertexContainer V, std::vector<SpannerCones> &H8, index_t start, index_t goal, number_t stretchBound ) {

            typedef std::pair<number_t,index_t>
                    DistanceIndexPair;
            typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
                    Heap;
            typedef Heap::handle_type
                    HeapHandle;

            size_t n = V.size();
            auto startPoint = V.at(start)->site().point();
            auto goalPoint = V.at(goal)->site().point();
            EuclideanDistanceToPoint h ( goalPoint ); // initialize heuristic functor

            Heap open;
            std::unordered_map<index_t,HeapHandle> handleToHeap(n);
            handleToHeap[start] = open.emplace( h( startPoint ), start );

            //unordered_set<size_t> closed(n);
            std::vector<index_t> parents(n);

            std::vector<number_t> g( n, INF );
            g[start] = 0;

            std::vector<number_t> f( n, INF );
            f[start] = h( startPoint );

            number_t baseEuclidean = f[start];

            DistanceIndexPair current = open.top(); // initialize current vertex to start
            index_t u_index = current.second;
            auto currentPoint = startPoint;
            auto neighborPoint = currentPoint;

            do {
                current = open.top();
                open.pop();

                u_index = current.second;
                currentPoint = V.at(u_index)->site().point();

                number_t currentStretch = f.at(u_index) / baseEuclidean;

                if (currentStretch > stretchBound)
                    return std::nullopt;


                if( u_index == goal ) return std::make_optional( g.at(goal) / baseEuclidean);

                for(auto cone : H8[u_index]) {

                    if (!cone.empty()) {

                        auto e = cone.front();

                        auto neighbor = (e.first == V.at(u_index)) ? e.second->storage_site().info() : e.first->storage_site().info();
                        neighborPoint = V.at(neighbor)->site().point();

                        if (u_index == start && neighbor == goal)
                            continue;

                        number_t newScore = g.at(u_index)
                                            + sqrt(CGAL::squared_distance( currentPoint, neighborPoint ));

                        if( newScore < g.at( neighbor ) ) {
                            parents[neighbor] = u_index;
                            g[neighbor] = newScore;
                            f[neighbor] = g.at(neighbor) + h(neighborPoint);
                            DistanceIndexPair q = std::make_pair( f.at(neighbor), neighbor );

                            if( contains( handleToHeap, neighbor ) ) {
                                HeapHandle neighborHandle = handleToHeap.at(neighbor);
                                open.update(neighborHandle,q);
                                open.update(neighborHandle);
                            } else {
                                handleToHeap[neighbor] = open.push(q);
                            }
                        }
                    }
                }

            } while( !open.empty() );

            return std::nullopt;
        }

        void decreaseDegree(std::vector<SpannerCones> &H8,
                            const std::vector<AnchorCones> &anchorEdges,
                            const std::vector<YaoCones> &yaoEdges,
                            const std::vector<FanCones> &pointFans,
                            const std::vector<NumYaoEdges> &yaoEdgeCount,
                            const std::vector<VertexHandle> &handles,
                            const DelaunayLinf &DT,
                            const VertexHandle u) {

            index_t u_id = u->storage_site().info();

            number_t stretchBound = MAX_STRETCH_CONSTRAINT;
            std::optional<SpannerEdge> candidateRemove = std::nullopt;

            for (cone_t cone = 0; cone < 4; cone++) {

                if (H8[u_id][cone].empty())
                    continue;

                SpannerEdge e = H8[u_id][cone].front();

                auto goal = (e.first == u) ? e.second->storage_site().info() : e.first->storage_site().info();

                auto stretch = AStar(handles, H8, u_id, goal, stretchBound);

                if (stretch) {

                    if (*stretch < stretchBound) {
                        stretchBound = *stretch;
                        candidateRemove = make_optional(e);
                    }
                }
            }

            // an edge is identified for removal
            if (candidateRemove) {

                auto p = (*candidateRemove).first;
                auto q = (*candidateRemove).second;

                index_t p_id = p->storage_site().info();
                index_t q_id = q->storage_site().info();

                SpannerEdge proposed = std::make_pair(p, q);
                SpannerEdge reversed = std::make_pair(q, p);


                for (auto &edges : H8[p_id]) {
                    edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
                    edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
                }

                for (auto &edges : H8[q_id]) {
                    edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
                    edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
                }

            }

            else {

                std::optional<SpannerEdge> candidateAdd = std::nullopt;
                stretchBound = MAX_STRETCH_CONSTRAINT;


                for (auto cone : H8[u_id]) {

                    auto e = cone.front();
                    index_t goal = (e.first == u) ? e.second->storage_site().info() : e.first->storage_site().info();

                    auto q = (e.first == u) ? e.second : e.first;
                    SpannerEdge edge = std::make_pair(u, q);

                    index_t q_id = q->storage_site().info();

                    VertexCirculator middle = DT.incident_vertices(q);

                    if (DT.is_infinite(middle)) {++middle;}
                    auto endpt = middle;

                    std::optional<number_t> stretch = std::nullopt;

                    do {

                        if (DT.is_infinite(middle)) {continue;}

                        index_t middle_id = middle->storage_site().info();
                        size_t currentCone = getCone(q, middle);

                        bool actualEdge = inEdgeList(H8[q_id][currentCone], q, middle) || inEdgeList(H8[middle_id][(currentCone+2)%4], q, middle);
                        size_t neighborDegree = accumulate(H8[middle_id].begin(), H8[middle_id].end(), 0, [](size_t sum, auto value) {
                            return sum + value.size(); });

                        bool possible = !actualEdge;

                        if (possible) {

                            SpannerEdge proposed = std::make_pair(q, middle);

                            if (middle != pointFans[q_id][currentCone].first && middle != pointFans[q_id][currentCone].second) {

                                auto previous = middle,
                                        next = middle;
                                --previous;
                                ++next;

                                index_t previous_id = previous->storage_site().info();
                                index_t next_id = next->storage_site().info();

                                size_t previousCone = getCone(previous, middle);
                                size_t nextCone = getCone(next, middle);

                                bool edgePair = (yaoEdges[previous_id][previousCone].first == middle && yaoEdges[middle_id][(previousCone+2)%4].first != previous
                                                 && anchorEdges[previous_id][previousCone].first != middle && anchorEdges[middle_id][(previousCone+2)%4].first != previous)
                                                && (yaoEdges[next_id][nextCone].first == middle && yaoEdges[middle_id][(nextCone+2)%4].first != next
                                                    && anchorEdges[next_id][nextCone].first != middle && anchorEdges[middle_id][(nextCone+2)%4].first != next);

                                possible &= !edgePair;

                            }

                            if (possible) {

                                H8[q_id][currentCone].push_back(proposed);
                                H8[middle_id][(currentCone+2)%4].push_back(proposed);

                                stretch = AStar(handles, H8, u_id, goal, stretchBound);

                                if (stretch) {

                                    if (*stretch < stretchBound) {
                                        stretchBound = *stretch;
                                        candidateRemove = make_optional(edge);
                                        candidateAdd = make_optional(proposed);
                                    }
                                }

                                SpannerEdge reversed = std::make_pair(proposed.second, proposed.first);

                                for (auto &edges : H8[q_id]) {
                                    edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
                                    edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
                                }

                                for (auto &edges : H8[middle_id]) {
                                    edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
                                    edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
                                }
                            }
                        }

                    } while (++middle != endpt && !stretch);

                }

                if (candidateRemove) {

                    auto p = (*candidateRemove).first;
                    auto q = (*candidateRemove).second;

                    index_t p_id = p->storage_site().info();
                    index_t q_id = q->storage_site().info();

                    SpannerEdge standardRemove = std::make_pair(p, q);
                    SpannerEdge reversedRemove = std::make_pair(q, p);

                    for (auto &edges : H8[p_id]) {
                        edges.erase(remove(edges.begin(), edges.end(), standardRemove), edges.end());
                        edges.erase(remove(edges.begin(), edges.end(), reversedRemove), edges.end());
                    }

                    for (auto &edges : H8[q_id]) {
                        edges.erase(remove(edges.begin(), edges.end(), standardRemove), edges.end());
                        edges.erase(remove(edges.begin(), edges.end(), reversedRemove), edges.end());
                    }

                    auto v = ((*candidateAdd).first == q) ? (*candidateAdd).second : (*candidateAdd).first;
                    index_t v_id = v->storage_site().info();

                    size_t addedCone = getCone(q, v);

                    H8[q_id][addedCone].push_back(*candidateAdd);
                    H8[v_id][(addedCone+2)%4].push_back(*candidateAdd);

                    size_t v_degree = accumulate(H8[v_id].begin(), H8[v_id].end(), 0, [](size_t sum, auto value) {
                        return sum + value.size(); });

                    while (v_degree > 3) {
                        decreaseDegree(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, DT, v);
                        v_degree--;
                    }
                }
            }
        }

        void degree3fromH4(std::vector<SpannerCones> &H8,
                           const std::vector<AnchorCones> &anchorEdges,
                           const std::vector<YaoCones> &yaoEdges, const std::vector<FanCones> &pointFans,
                           const std::vector<NumYaoEdges> &yaoEdgeCount, const std::vector<VertexHandle> &handles,
                           const DelaunayLinf &DT) {

            for (const auto &u : handles) {

                index_t u_id = u->storage_site().info();
                size_t degree = accumulate(H8[u_id].begin(), H8[u_id].end(), 0, [](size_t sum, auto value) {
                    return sum + value.size(); });

                if (degree > 3) {

                    decreaseDegree(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, DT, u);

                }
            }
//
//                    number_t stretchBound = MAX_STRETCH_CONSTRAINT;
//                    optional<SpannerEdge> candidateRemove = nullopt;
//
//                    for (cone_t cone = 0; cone < 4; cone++) {
//
//                        if (H8[u_id][cone].size() == 0)
//                            continue;
//
//                        SpannerEdge e = H8[u_id][cone].front();
//
//                        auto goal = (e.first == u) ? e.second->storage_site().info() : e.first->storage_site().info();
//
//                        auto stretch = AStar(handles, H8, u_id, goal, stretchBound);
//
//                        if (stretch) {
//
//                            if (*stretch < stretchBound) {
//                                stretchBound = *stretch;
//                                candidateRemove = make_optional(e);
//                            }
//                        }
//                    }
//
//                    // an edge is identified for removal
//                    if (candidateRemove) {
//
//                        auto p = (*candidateRemove).first;
//                        auto q = (*candidateRemove).second;
//
//                        index_t p_id = p->storage_site().info();
//                        index_t q_id = q->storage_site().info();
//
//                        SpannerEdge proposed = std::make_pair(p, q);
//                        SpannerEdge reversed = std::make_pair(q, p);
//
//
//                        for (auto &edges : H8[p_id]) {
//                            edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
//                            edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
//                        }
//
//                        for (auto &edges : H8[q_id]) {
//                            edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
//                            edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
//                        }
//
//                    }
//
//                    else {
//
//                        optional<SpannerEdge> candidateAdd = nullopt;
//                        number_t stretchBound = MAX_STRETCH_CONSTRAINT;
//
//                        for (auto cone : H8[u_id]) {
//
//                            auto e = cone.front();
//                            index_t goal = (e.first == u) ? e.second->storage_site().info() : e.first->storage_site().info();
//
//                            auto q = (e.first == u) ? e.second : e.first;
//                            SpannerEdge edge = std::make_pair(u, q);
//
//                            index_t q_id = q->storage_site().info();
//
//                            VertexCirculator middle = DT.incident_vertices(q);
//
//                            if (DT.is_infinite(middle)) {++middle;}
//                            auto endpt = middle;
//
//                            optional<number_t> stretch = nullopt;
//
//                            do {
//
//                                if (DT.is_infinite(middle)) {continue;}
//
//                                index_t middle_id = middle->storage_site().info();
//                                size_t currentCone = getCone(q, middle);
//
//                                bool actualEdge = inEdgeList(H8[q_id][currentCone], q, middle) || inEdgeList(H8[middle_id][(currentCone+2)%4], q, middle);
//                                size_t neighborDegree = accumulate(H8[middle_id].begin(), H8[middle_id].end(), 0, [](size_t sum, auto value) {
//                                    return sum + value.size(); });
//
//                                bool possible = !actualEdge && neighborDegree < 3;
//
//                                if (possible) {
//
//                                    SpannerEdge proposed = std::make_pair(q, middle);
//
//                                    if (middle != pointFans[q_id][currentCone].first && middle != pointFans[q_id][currentCone].second) {
//
//                                        auto previous = middle,
//                                                next = middle;
//                                        --previous;
//                                        ++next;
//
//                                        index_t previous_id = previous->storage_site().info();
//                                        index_t next_id = next->storage_site().info();
//
//                                        size_t previousCone = getCone(previous, middle);
//                                        size_t nextCone = getCone(next, middle);
//
//                                        bool edgePair = (yaoEdges[previous_id][previousCone].first == middle && yaoEdges[middle_id][(previousCone+2)%4].first != previous
//                                                         && anchorEdges[previous_id][previousCone].first != middle && anchorEdges[middle_id][(previousCone+2)%4].first != previous)
//                                                        && (yaoEdges[next_id][nextCone].first == middle && yaoEdges[middle_id][(nextCone+2)%4].first != next
//                                                            && anchorEdges[next_id][nextCone].first != middle && anchorEdges[middle_id][(nextCone+2)%4].first != next);
//
//                                        possible &= !edgePair;
//
//                                    }
//
//                                    if (possible) {
//
//                                        H8[q_id][currentCone].push_back(proposed);
//                                        H8[middle_id][(currentCone+2)%4].push_back(proposed);
//
//                                        stretch = AStar(handles, H8, u_id, goal, stretchBound);
//
//                                        if (stretch) {
//
//                                            if (*stretch < stretchBound) {
//                                                stretchBound = *stretch;
//                                                candidateRemove = make_optional(edge);
//                                                candidateAdd = make_optional(proposed);
//                                            }
//                                        }
//
//                                        SpannerEdge reversed = std::make_pair(proposed.second, proposed.first);
//
//                                        for (auto &edges : H8[q_id]) {
//                                            edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
//                                            edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
//                                        }
//
//                                        for (auto &edges : H8[middle_id]) {
//                                            edges.erase(remove(edges.begin(), edges.end(), proposed), edges.end());
//                                            edges.erase(remove(edges.begin(), edges.end(), reversed), edges.end());
//                                        }
//                                    }
//                                }
//
//                            } while (++middle != endpt && !stretch);
//
//                        }
//
//                        if (candidateRemove) {
//
//                            auto p = (*candidateRemove).first;
//                            auto q = (*candidateRemove).second;
//
//                            index_t p_id = p->storage_site().info();
//                            index_t q_id = q->storage_site().info();
//
//                            SpannerEdge standardRemove = std::make_pair(p, q);
//                            SpannerEdge reversedRemove = std::make_pair(q, p);
//
//                            for (auto &edges : H8[p_id]) {
//                                edges.erase(remove(edges.begin(), edges.end(), standardRemove), edges.end());
//                                edges.erase(remove(edges.begin(), edges.end(), reversedRemove), edges.end());
//                            }
//
//                            for (auto &edges : H8[q_id]) {
//                                edges.erase(remove(edges.begin(), edges.end(), standardRemove), edges.end());
//                                edges.erase(remove(edges.begin(), edges.end(), reversedRemove), edges.end());
//                            }
//
//                            auto v = ((*candidateAdd).first == q) ? (*candidateAdd).second : (*candidateAdd).first;
//                            index_t v_id = v->storage_site().info();
//
//                            size_t addedCone = getCone(q, v);
//
//                            H8[q_id][addedCone].push_back(*candidateAdd);
//                            H8[v_id][(addedCone+2)%4].push_back(*candidateAdd);
//
//                        }
//                    }
//
//                }

            // the greedy addition of the edges after edge removal
            for (auto u : handles) {

                index_t u_id = u->storage_site().info();
                size_t degree = accumulate(H8[u_id].begin(), H8[u_id].end(), 0, [](size_t sum, auto value) {
                    return sum + value.size(); });

                bool found = true;
                bool addEdges = degree < 3;

                while (addEdges) {

                    found = false;

                    VertexCirculator v = DT.incident_vertices(u);
                    if (DT.is_infinite(v))
                        ++v;

                    index_t v_id = v->storage_site().info();
                    number_t min_distance = INF;

                    std::optional<VertexHandle> optimal = std::nullopt;

                    auto endpt = v;

                    do {

                        if (DT.is_infinite(v))
                            continue;

                        v_id = v->storage_site().info();
                        size_t v_degree = accumulate(H8[v_id].begin(), H8[v_id].end(), 0, [](size_t sum, auto value) {
                            return sum + value.size(); });

                        size_t currentCone = getCone(u,v);

                        bool alreadyEdge = inEdgeList(H8[u_id][currentCone], u, v) || inEdgeList(H8[v_id][(currentCone+2)%4], u, v);
                        bool possible = !alreadyEdge && (v_degree < 3);

                        if (possible) {

                            bool shortcut = false;
                            auto fan = pointFans[u_id][currentCone];

                            if (v != fan.first && v != fan.second) {

                                auto previous = v,
                                        middle = v,
                                        next = v;

                                --previous;
                                ++next;

                                index_t previous_id = previous->storage_site().info();
                                index_t next_id = next->storage_site().info();
                                index_t middle_id = middle->storage_site().info();

                                size_t previousCone = getCone(previous, middle);
                                size_t nextCone = getCone(next, middle);

                                shortcut  = (yaoEdges[previous_id][previousCone].first == middle && yaoEdges[middle_id][(previousCone+2)%4].first != previous
                                                  && anchorEdges[previous_id][previousCone].first != middle && anchorEdges[middle_id][(previousCone+2)%4].first != previous)
                                                 && (yaoEdges[next_id][nextCone].first == middle && yaoEdges[middle_id][(nextCone+2)%4].first != next
                                                     && anchorEdges[next_id][nextCone].first != middle && anchorEdges[middle_id][(nextCone+2)%4].first != next);

                            }

                            if (!shortcut) {

                                CGAL::Point_2 u_point = u->site().point();
                                CGAL::Point_2 v_point = v->site().point();

                                auto v_distance = (number_t) sqrt(CGAL::squared_distance(u_point, v_point));

                                if (v_distance < min_distance) {

                                    v_id = v->storage_site().info();

                                    found = true;
                                    min_distance = v_distance;
                                    optimal = std::make_optional(v);

                                }
                            }
                        }

                    } while (++v != endpt);

                    if (found) {

                        cone_t newCone = getCone(u, *optimal);
                        SpannerEdge newEdge = std::make_pair(u, *optimal);

                        index_t optimal_id = (*optimal)->storage_site().info();

                        H8[u_id][newCone].push_back(newEdge);
                        H8[optimal_id][(newCone+2)%4].push_back(newEdge);

                    }

                    degree = accumulate(H8[u_id].begin(), H8[u_id].end(), 0, [](size_t sum, auto value) {
                        return sum + value.size(); });

                    addEdges = (degree < 3) && found;
                }
            }
        }
    }





    void DEG3(const bdps::input_t& in, bdps::output_t& out) {

        using namespace degree3;
        using degree3::VertexHandle, degree3::VertexCirculator, degree3::FaceHandle;

        const index_t n = in.size();
        if (n > SIZE_T_MAX - 1 || n <= 1) return;

        // construct Linf Delaunay triangulation
        std::vector<Point> P(in);
        std::vector<size_t> index;
        spatialSort<K>(P, index);
        DelaunayLinf DT;
        DelaunayLinf::Site_2 site;
        index_t id = 0;

        // store the vertex handles
        std::vector<VertexHandle> handles(n);

        //FaceHandle hint;
        for (size_t entry : index) {
            Point p = P[entry];
            site = DelaunayLinf::Site_2::construct_site_2(p);
            auto vh = DT.insert(site, entry);
            //hint = vh->face();
            //vh->storage_site().info() = entry;
            handles[entry] = vh;
        }
        //construct YaoEdges
        std::vector<YaoCones> yaoEdges(n, YaoCones(4));
        std::vector<FanCones> pointFans(n, FanCones(4));
        std::vector<NumYaoEdges> yaoEdgeCount(n, NumYaoEdges(4));
        addYaoEdges(yaoEdges, pointFans, yaoEdgeCount, handles, DT);

        // identify the anchors
        std::vector<AnchorCones> anchorEdges(n, AnchorCones(4, std::make_pair(DT.infinite_vertex(), None)));
        determineAnchors(anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, DT);

        // construct H8 --> degree 8 spanner
        std::vector<SpannerCones> H8(n, SpannerCones(4));
        degreeEightSpanner(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, DT);

        // prune edges from H8 to produce H4 --> degree 4 spanner
        processSpanner(H8, anchorEdges, yaoEdges, pointFans, handles, DT);

        // produce degree 3 spanner from H4
        degree3fromH4(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, DT);

        std::vector<index_tPair> edgeList;

        for (const auto &u : handles) {
            for (size_t cone = 0; cone < 4; cone++) {
                for (const auto &edge : H8[u->storage_site().info()][cone]) {
                    edgeList.emplace_back((edge.first)->storage_site().info(),
                                          (edge.second)->storage_site().info());
                }
            }
        }

        std::copy(edgeList.begin(), edgeList.end(), std::back_inserter(out));
    }


}



#endif //SPANNERS_DEGREE3_H