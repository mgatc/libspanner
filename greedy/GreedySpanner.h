//
// Created by justin on 4/2/22.
//

#ifndef SP_MER_17_GREEDYSPANNER_H
#define SP_MER_17_GREEDYSPANNER_H

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
    double shortestPath(const std::vector<Point> &points,
                        std::pair<int, int> cur,
                        const std::unordered_map<int, std::vector<int>> &adjMap) {
        std::vector<number_t> distances(points.size(), DBL_MAX);
        distances[cur.first] = 0.0;
        //priority_queue<pair<number_t, int>, vector<pair<number_t,int>>, greater<>> pq;
        std::map<number_t, int> lookup;
        std::vector<std::map<number_t, int>::iterator> handles(points.size(), lookup.end());
        //pq.push(make_pair(0.0, cur.first));
        auto x = lookup.emplace(0.0, cur.first);
        handles[cur.first] = x.first;
        while (!lookup.empty()) {
            auto curr = lookup.begin();
            auto adj = adjMap.at(curr->second);
            for (auto a: adj) {
                double pathTotal_With_Adjacent_Distance = curr->first + getDistance(points[curr->second], points[a]);
                //distance_to_current + getDistance(points[current_position], points[a]);
                if (distances[a] > pathTotal_With_Adjacent_Distance) {
                    distances[a] = pathTotal_With_Adjacent_Distance;
                    //pq.push(make_pair(pathTotal_With_Adjacent_Distance, a));
                    if (handles[a] != lookup.end()) lookup.erase(handles[a]);
                    auto y = lookup.emplace(pathTotal_With_Adjacent_Distance, a);
                    handles[a] = y.first;
                }
            }
            lookup.erase(curr);
        }
        return distances[cur.second];
    }

    void sortPairs(const std::vector<Point> &points,
                   std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::function<bool(
                           std::pair<int, int>, std::pair<int, int>)>> &pq) {
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                pq.push(std::make_pair(i, j));
            }
        }
    }
}

    void GreedySpanner(const greedy::input_t& points, greedy::output_t& out, const double t = 1.0) {
//    void GreedySpanner(const std::vector<Point> &points, const std::vector<int> &indices,
//                       std::vector<Edge> &edges, double t,
//                       std::unordered_map<int, std::vector<int>> &adjMap){
        assert( t > 1.0 - std::numeric_limits<double>::epsilon());

        std::function<bool(std::pair<int,int>, std::pair<int,int>)> cmp = [&points](std::pair<int,int> uv, std::pair<int,int> wx){
            return (getDistance(points[uv.first], points[uv.second])) > (getDistance(points[wx.first], points[wx.second]));
        };
        std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> pairs(cmp);

        classic_greedy::sortPairs(points, pairs);
        std::unordered_map<int, std::vector<int>> adjMap;
        while(!pairs.empty()){
            if(classic_greedy::shortestPath(points, pairs.top(), adjMap) >
               t * getDistance(points[pairs.top().first], points[pairs.top().second])){
                out.emplace_back(std::make_pair(pairs.top().first, pairs.top().second));
                adjMap[pairs.top().first].emplace_back(pairs.top().second);
                adjMap[pairs.top().second].emplace_back(pairs.top().first);
            }
            pairs.pop();
        }
    }


}

#endif //SP_MER_17_GREEDYSPANNER_H
