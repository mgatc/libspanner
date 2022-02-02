#ifndef UNF_SPANNERS_BKPX2015_H
#define UNF_SPANNERS_BKPX2015_H

#include <array>
#include <iostream>
#include <list>
#include <cassert>
#include <string>
#include <limits>
#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles
#include <boost/functional/hash.hpp> // size_t pair hash : used in Yao_inf_4

#include <CGAL/algorithm.h> //
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/circulator.h>

#include "constants.h"
#include "delaunay/DelaunayLinf.h"
#include "../bdps/types.h"
#include "Utilities.h"

namespace spanner {


    namespace bkpx2015 {

        typedef DelaunayLinf::Site_2 Site;
        //typedef Site::Point_2 Point;
        typedef DelaunayLinf::Vertex_circulator VertexCirculator;
        //typedef DelaunayLinf::Edge LinfEdge;
        typedef DelaunayLinf::Vertex_handle VertexHandle;
        typedef DelaunayLinf::Face_handle FaceHandle;

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

                        number_t proposedDistance = CGAL::l_infinity_distance(point->site().point(),
                                                                              circ->site().point());

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

                } while (++circ !=
                         endpoint); // finished determining the Yao edges + how many points are in fan of u's cone i
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
                                            (yaoEdges[vhigher->storage_site().info()][getCone(vhigher, v)].first ==
                                             v) &&
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
                                                       yaoEdges[current->storage_site().info()][getCone(current,
                                                                                                        previous)].first ==
                                                       previous &&
                                                       yaoEdges[previous->storage_site().info()][getCone(previous,
                                                                                                         current)].first !=
                                                       current);

                                inCanonical = inFan && unidirectional;

                                bool yaoConnected = inCanonical && (yaoEdges[u_id][cone].first == current ||
                                                                    yaoEdges[current->storage_site().info()][
                                                                            (cone + 2) % 4].first == u);

                                if (yaoConnected)
                                    crown = current;

                            }

                            assert(yaoEdges[u_id][cone].first == crown ||
                                   yaoEdges[crown->storage_site().info()][(cone + 2) % 4].first == u);

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
                                               || DT.is_infinite(
                            anchorEdges[v->storage_site().info()][(cone + 2) % 4].first))) {
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

                            inChain = !(anchorEdges[previous_id][currentCone].second == Strong ||
                                        anchorEdges[previous_id][currentCone].second == StrongSelected);

                            if (!inChain)
                                break;

                            current = anchorEdges[previous_id][currentCone].first;

                            assert(!DT.is_infinite(current));

                        } while (inChain);

                        assert(anchorEdges[previous_id][currentCone].second == Strong ||
                               anchorEdges[previous_id][currentCone].second == StrongSelected);

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
    }


    void BKPX2015(const bdps::input_t &in, bdps::output_t &out) {

        using namespace bkpx2015;

        using bkpx2015::VertexHandle, bkpx2015::VertexCirculator, bkpx2015::FaceHandle;

        // construct Linf Delaunay triangulation
        std::vector<Point> P(in);
        std::vector<size_t> index;
        spatialSort<K>(P, index);
        DelaunayLinf DT;
        Site site;
        index_t id = 0;

        // store the vertex handles
        const size_t n = P.size();
        std::vector<VertexHandle> handles(n);

        //FaceHandle hint;
        for (size_t entry : index) {
            Point p = P[entry];
            site = Site::construct_site_2(p);
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

        processSpanner(H8, anchorEdges, yaoEdges, pointFans, handles, DT);

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
        // Send resultant graph to output iterator
//        for (auto e : edgeList) {
//            // Edge list is only needed for printing. Remove for production.
//            //edgeList.emplace_back(handles.at(e.first)->point(), handles.at(e.second)->point());
//
//            *result = e;
//            ++result;
//        }

        // START PRINTER NONSENSE
//    if(printLog) {
//        GraphPrinter printer(1.25); // argument number is scaling factor --> manipulate based on size of point set
//        GraphPrinter::OptionsList options;
//
//        options = {
//            {"color", printer.inactiveEdgeColor},
//            {"line width", to_string(printer.inactiveEdgeWidth)}
//        };
////        printer.drawEdgesOfSDG(DT, options);
//
//        options = { // active edge options
//            {"color", printer.activeEdgeColor},
//            {"line width", to_string(printer.activeEdgeWidth)}
//        };
//
//        vector<pair<Point, Point>> pointEdgeList;
//        pointEdgeList.reserve(edgeList.size());
//
//        for (auto e : edgeList) {
//
//            pointEdgeList.emplace_back(P.at(e.first), P.at(e.second));
//
//        }
//
//        printer.drawEdges(edgeList.begin(), edgeList.end(), P, options);
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
//        printer.drawVerticesWithInfoSDG(DT, options, borderOptions);
//
//        printer.print("BKPX2015");
//        cout << "\n";
//    }
        // END PRINTER NONSENSE


    }


}


#endif // UNF_SPANNERS_BKPX2015_H
