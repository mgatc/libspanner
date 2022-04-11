//
// Created by justin on 4/2/22.
//

#ifndef SP_MER_17_GREEDYSPANNER_H
#define SP_MER_17_GREEDYSPANNER_H

#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>


#include "../geometry.h"
#include "../greedy/types.h"
#include "../types.h"
#include "../utilities.h"

namespace spanner {
namespace classic_greedy {
    number_t shortestPath(const std::vector<Point> &P,
                        const Edge &currentPair,
                        const AdjacencyListDense &adjMap) {

        std::unordered_map<index_t,number_t> g;
        std::unordered_map<index_t,number_t> f;

        g[currentPair.first] = 0.0;
        f[currentPair.first] = getDistance(P.at(currentPair.first), P.at(currentPair.second));

        std::map<number_t, index_t> open;
        std::unordered_map<index_t,std::map<number_t, index_t>::iterator> handles;

        auto x = open.emplace(f.at(currentPair.first), currentPair.first);
        handles[currentPair.first] = x.first;

        while (!open.empty()) {
            auto curr = open.begin();
            if(curr->second == currentPair.second) {
                return g.at(currentPair.second);
            }
//            auto& adj = adjMap.at(curr->second);
            for (auto a : adjMap.at(curr->second)) {
                number_t newG = g.at(curr->second) + getDistance(P.at(curr->second), P.at(a));
                //distance_to_current + getDistance(P[current_position], P[a]);
                if (!contains(g, a) || g.at(a) > newG) {
                    g[a] = newG;
                    f[a] = newG + getDistance(P[a], P[currentPair.second]);
                    //pq.push(make_pair(newG, a));
                    if (contains(handles,a) ) open.erase(handles[a]);
                    auto y = open.emplace(f.at(a), a);
                    handles[a] = y.first;
                }
            }
            open.erase(curr);
        }
        return INF;
    }

}

    void GreedySpanner(const greedy::input_t& P, greedy::output_t& out, const double t = 1.0) {
//    void GreedySpanner(const std::vector<Point> &P, const std::vector<int> &indices,
//                       std::vector<Edge> &edges, double t,
//                       std::unordered_map<int, std::vector<int>> &adjMap){
        assert( t > 1.0 - std::numeric_limits<double>::epsilon());

        const index_t n = P.size();

        AdjacencyListDense adjMap(n, std::unordered_set<index_t>());

        std::vector<Edge> allPairs;
        allPairs.reserve(n*n);

        for (size_t i = 0; i < P.size(); ++i) {
            for (size_t j = i + 1; j < P.size(); ++j) {
                allPairs.emplace_back(i, j);
            }
        }
        std::sort(allPairs.begin(), allPairs.end(), [&P]( const Edge& uv,  const Edge& wx){
            return (CGAL::squared_distance(P[uv.first], P[uv.second])) < (CGAL::squared_distance(P[wx.first], P[wx.second]));
        });

        for(auto currentPair : allPairs ){
            if(classic_greedy::shortestPath(P, currentPair, adjMap) >
               static_cast<number_t>(t) * getDistance(P[currentPair.first], P[currentPair.second])){
                out.emplace_back(currentPair.first, currentPair.second);
                adjMap[currentPair.first].insert(currentPair.second);
                adjMap[currentPair.second].insert(currentPair.first);
            }
        }
    }


}

#endif //SP_MER_17_GREEDYSPANNER_H
