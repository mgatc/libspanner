#ifndef LIBSPANNER_DEGREE_H
#define LIBSPANNER_DEGREE_H

#include <algorithm> // swap, max_element, accumulate
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../utilities.h"

namespace spanner {

    template<typename RandomAccessIterator>
    size_t degree(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        const std::vector<EdgeType> edges(edgesBegin, edgesEnd);

        if(edges.empty()) {
            return 0;
        }

        std::unordered_map<VertexType, std::unordered_set<VertexType>> adj;
        // for each edge
        for (auto e : edges) {
            auto first = adj.begin();
            std::tie(first, std::ignore) = adj.emplace(e.first, std::unordered_set<VertexType>());
            (*first).second.insert(e.second);

            auto second = adj.begin();
            std::tie(second, std::ignore) = adj.emplace(e.second, std::unordered_set<VertexType>());
            (*second).second.insert(e.first);
        }
        auto max_el = std::max_element(adj.begin(), adj.end(), [&](const auto &lhs, const auto &rhs) {
            return lhs.second.size() < rhs.second.size();
        });
//        cout<<"Largest degree vertex="<<max_el->first<<endl;

        return max_el->second.size();
    }

    template<typename RandomAccessIterator>
    number_t avgDegreePerPoint(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        const std::vector<EdgeType> edges(edgesBegin, edgesEnd);
        if(edges.empty()) {
            return 0;
        }
        std::unordered_map<VertexType, std::unordered_set<VertexType>> adj;
        // for each edge
        for (auto e : edges) {
            auto first = adj.begin();
            std::tie(first, std::ignore) = adj.emplace(e.first, std::unordered_set<VertexType>());
            (*first).second.insert(e.second);

            auto second = adj.begin();
            std::tie(second, std::ignore) = adj.emplace(e.second, std::unordered_set<VertexType>());
            (*second).second.insert(e.first);
        }
        auto avg = std::accumulate(adj.begin(), adj.end(), 0.0, [&](const number_t &sum, const auto &current) {
            return sum + current.second.size();
        }) / number_t(adj.size());

        return avg;
    }

//    template<typename Triangulation>
//    size_t degree(const Triangulation &T) {
//        typedef typename Triangulation::Point
//                Point_2;
//
//        // fill a vector with edges so we can call the range-based degree function
//        const std::vector<std::pair<Point_2, Point_2>> edges;
//        edges.reserve(T.number_of_vertices());
//
//        for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
//            auto p = make_pair(
//                    e->first->vertex((e->second + 1) % 3)->point(),
//                    e->first->vertex((e->second + 2) % 3)->point()
//            );
//            // Add both in and out edges
//            forBoth(p, [&](Point a, Point b) {
//                edges.emplace_back(a, b);
//            });
//        }
//        return degree(edges.begin(), edges.end());
//    }

} // namespace spanner


#endif // LIBSPANNER_DEGREE_H


