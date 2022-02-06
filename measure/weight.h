#ifndef LIBSPANNER_WEIGHT_H
#define LIBSPANNER_WEIGHT_H

#include <tuple>

#include "../geometry.h"
#include "../points/mst.h"
#include "../types.h"
#include "../utilities.h"

namespace spanner {


    template< class VertexIterator, class EdgeIterator>
    number_t weight( VertexIterator pointsBegin,
                     VertexIterator pointsEnd,
                     EdgeIterator edgesBegin,
                     EdgeIterator edgesEnd ) {
        const std::vector<Point> P(pointsBegin,pointsEnd);

        number_t w = 0.0;
        index_t p,q;
        for (auto e = edgesBegin; e != edgesEnd; ++e) {
            std::tie(p,q) = *e;
            w += getDistance(P[p],P[q]);
        }
        return w;
    }

//    template<typename Triangulation>
//    number_t weight(const Triangulation &T) {
//        number_t w = 0;
//        for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
//            auto p = std::make_pair(
//                    e->first->vertex((e->second + 1) % 3)->point(),
//                    e->first->vertex((e->second + 2) % 3)->point()
//            );
//            w += getDistance(p.first, p.second);
//        }
//        return w;
//    }

    template< class VertexIterator, class EdgeIterator>
    number_t getLightness( VertexIterator pointsBegin,
                           VertexIterator pointsEnd,
                           EdgeIterator edgesBegin,
                           EdgeIterator edgesEnd ) {
        std::vector<Point> P(pointsBegin,pointsEnd);
        std::vector<Edge> E(edgesBegin,edgesEnd);
        std::list<Edge> MST;
        getMST( P.begin(), P.end(), E.begin(), E.end(), back_inserter(MST) );
        number_t weightOfMST = weight(P.begin(), P.end(), MST.begin(), MST.end() ),
                weightOfG   = weight(P.begin(), P.end(), E.begin(), E.end() ),
                lightness = weightOfG / weightOfMST;
        return lightness;
    }

} // namespace spanner

#endif // LIBSPANNER_WEIGHT_H


