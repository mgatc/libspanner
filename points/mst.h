#ifndef LIBSPANNER_MST_H
#define LIBSPANNER_MST_H

#include <list>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include "../geometry.h"
#include "../types.h"

namespace spanner {

    template< class VertexIterator, class EdgeIterator, class EdgeOutputIterator >
    void getMST( VertexIterator pointsBegin,
                 VertexIterator pointsEnd,
                 EdgeIterator edgesBegin,
                 EdgeIterator edgesEnd,
                 EdgeOutputIterator out)
    {

        using namespace boost;
        typedef adjacency_list<vecS, vecS, undirectedS,
                Point,
                property<edge_weight_t,number_t>
        > Graph;

        Graph G;

        for( auto pit=pointsBegin; pit!=pointsEnd; ++pit ) {
            auto v = add_vertex( *pit, G );
        }

        index_t p, q;
        for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
            std::tie( p, q ) = *eit;
            number_t wt = getDistance( G[p], G[q] );
            auto e = add_edge( p, q, wt, G );
        }

        typedef Graph::vertex_descriptor VertexDescriptor;
        typedef Graph::edge_descriptor EdgeDescriptor;

        std::list<EdgeDescriptor> mst;
        boost::kruskal_minimum_spanning_tree(G,std::back_inserter(mst));

        for(auto it = mst.begin(); it != mst.end(); ++it){
            EdgeDescriptor ed = *it;
            VertexDescriptor u = source(ed, G),
                    v = target(ed, G);
            *out = std::make_pair(u,v);
        }
    }

} // namespace spanner


#endif // LIBSPANNER_MST_H


