//
// Created by matt on 2/6/22.
//

#ifndef LIBSPANNER_MEASURE_H
#define LIBSPANNER_MEASURE_H

#include <iostream>
#include <list>
#include <utility>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/ch_selected_extreme_points_2.h>

#include "../delaunay/DelaunayL2.h"
#include "../geometry.h"
#include "../types.h"

namespace spanner {

    std::pair<Point,Point>
    getBoundingBox(const PointContainer &in) {
        using std::make_pair;
        using std::cout;

        if( in.empty() )
            return make_pair( Point(0,0), Point(0,0) );

        auto x_min = in.begin(),
                x_max = x_min,
                y_min = x_min,
                y_max = x_min;

        CGAL::ch_nswe_point( in.begin(), in.end(),
                             y_max,
                             y_min,
                             x_min,
                             x_max );
        auto bbox = make_pair( Point(x_min->x(),y_min->y()),
                               Point(x_max->x(),y_max->y()) );
        //cout<<" bbox: "<<bbox.first<< " " <<bbox.second<<"\n";
        return bbox;
    }
    number_t getBoundingBoxLength(const PointContainer &in) {
        auto bbox = getBoundingBox(in);
        return abs( bbox.first.y() - bbox.second.y() );
    }
    number_t getBoundingBoxWidth(const PointContainer &in) {
        auto bbox = getBoundingBox(in);
        return abs( bbox.first.x() - bbox.second.x() );
    }
    number_t getBoundingBoxArea(const PointContainer &in) {
        return getBoundingBoxLength(in) * getBoundingBoxWidth(in);
    }

    size_t getConvexHullSize(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        if(T.number_of_vertices() <= 3)
            return T.number_of_vertices();
        using std::cout;
        auto hull = T.incident_vertices( T.infinite_vertex() ),
                done = hull;
        size_t size = 0;
        //cout<<" walking hull...\n";
        do {
            ++size;
            //    cout<<"  "<<size<<") "<<hull->point()<<"\n";
        } while( ++hull != done );

        return size;
    }
    number_t getConvexHullWeight(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        if(T.number_of_vertices() < 2)
            return 0.0;
        using std::cout;
        auto hull = T.incident_vertices( T.infinite_vertex() ),
                prev = hull,
                done = hull;
        --prev;

        number_t weight = 0;
        size_t size = 0;
        //cout<<" walking hull...\n";
        do {
            number_t edgeWt = CGAL::sqrt( CGAL::squared_distance(prev->point(), hull->point()) );
            weight += edgeWt;
            ++size;

            //cout<< "  " << size << ") " << hull->point() <<" "<< edgeWt<< " " << weight << "\n";
        } while( ++hull != done );

        return weight;
    }
    number_t closestPair(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        using std::cout;

        //std::pair<size_t,size_t> closestPair = std::make_pair(0,1);
        number_t closestDistanceSquared = INF,
                currentDistanceSquared = INF;
        //cout<<" finding closest pair...\n";
        for( auto eit=T.finite_edges_begin();
             eit!=T.finite_edges_end(); ++eit ) {
            auto edge = T.segment(*eit);
            currentDistanceSquared = edge.squared_length();

            if( currentDistanceSquared < closestDistanceSquared )
                closestDistanceSquared = currentDistanceSquared;
        }

        number_t closestDistance = CGAL::sqrt( closestDistanceSquared );

        return closestDistance;
    }
    number_t farthestPair(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        if(T.number_of_vertices()<2)
            return 0.0;
        using std::cout;
        std::vector<Point> P;

        auto hull = T.incident_vertices( T.infinite_vertex() ),
                done = hull;

        do {
            P.push_back( hull->point() );
        } while( ++hull != done );

        number_t d_max_sq = 0;

        const size_t n = P.size();

        //cout<<" finding farthest pair...\n";

        for( size_t p=0; p<n; ++p ) {
            for( size_t q=p+1; q<n; ++q ) {
                number_t current = CGAL::squared_distance(P[p],P[q]);
                if( current > d_max_sq ) {
                    d_max_sq = current;
                    //cout<<"  "<<p<<" "<<q<<" "<<CGAL::sqrt(current)<<"\n";
                }
            }
        }
        return CGAL::sqrt( d_max_sq );
    }

    number_t spread(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        if(T.number_of_vertices() < 2)
            return 0.0;
        using std::cout;
        //cout<<"  calculating spread...";
        number_t farthest = farthestPair(in),
                closest = closestPair(in),
                spread = farthest / closest;

        //cout<<spread<<"\n";

        return spread;
    }

    template <typename T>
    struct Is_finite {
        const T* t_;
        Is_finite()
                : t_(NULL)
        {}
        Is_finite(const T& t)
                : t_(&t)
        { }
        template <typename VertexOrEdge>
        bool operator()(const VertexOrEdge& voe) const {
            return ! t_->is_infinite(voe);
        }
    };

    number_t weightOfMST(const PointContainer &in) {
        DelaunayL2 T(in.begin(),in.end());
        if(T.number_of_vertices() < 2)
            return 0.0;

        using std::cout;
        using namespace boost;
        //cout<<" determining weight of MST...";

        typedef Is_finite<DelaunayL2> Filter;
        typedef boost::filtered_graph<DelaunayL2,Filter,Filter> Finite_triangulation;
        typedef typename boost::graph_traits<Finite_triangulation>::vertex_descriptor vertex_descriptor;
        typedef typename boost::graph_traits<Finite_triangulation>::vertex_iterator vertex_iterator;
        typedef typename boost::graph_traits<Finite_triangulation>::edge_descriptor edge_descriptor;
        // The BGL makes use of indices associated to the vertices
        // We use a std::map to store the index
        typedef std::map<vertex_descriptor,int> VertexIndexMap;
        // A std::map is not a property map, because it is not lightweight
        typedef boost::associative_property_map<VertexIndexMap> VertexIndexPMap;
        // Associate indices to the vertices
        VertexIndexMap vertex_id_map;
        VertexIndexPMap vertex_index_pmap(vertex_id_map);

        Filter is_finite(T);
        Finite_triangulation ft(T, is_finite, is_finite);

        vertex_iterator vit, ve;
        // Associate indices to the vertices
        int index = 0;
        // boost::tie assigns the first and second element of the std::pair
        // returned by boost::vertices to the variables vit and ve
        for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
            vertex_descriptor  vd = *vit;
            vertex_id_map[vd]= index++;
        }
        // We use the default edge weight which is the squared getBoundingBoxLength of the edge
        // This property map is defined in graph_traits_Triangulation_2.h
        // In the function call you can see a named parameter: vertex_index_map
        std::list<edge_descriptor> mst;
        boost::kruskal_minimum_spanning_tree(ft,
                                             std::back_inserter(mst),
                                             vertex_index_map(vertex_index_pmap));

        number_t wt = 0.0;
        for(auto it = mst.begin(); it != mst.end(); ++it){
            edge_descriptor ed = *it;
            vertex_descriptor svd = source(ed,T);
            vertex_descriptor tvd = target(ed,T);
            DelaunayL2::Vertex_handle sv = svd;
            DelaunayL2::Vertex_handle tv = tvd;
            number_t edgeWt = CGAL::sqrt( CGAL::squared_distance( sv->point(), tv->point() ) );
            wt += edgeWt;
            //std::cout << "  " << sv->point() << " | " << tv->point() << " " << edgeWt<<" "<< weight<<"\n";
        }
        //cout<<weight<<"\n";
        return wt;
    }

    Point midpoint(const PointContainer &in) {
        if(in.empty())
            return Point(0.0,0.0);

        DelaunayL2 T(in.begin(),in.end());
        number_t x = 0,
                y = 0;

        for( auto p : in )
        {
            x += p.x();
            y += p.y();
        }
        x /= static_cast<number_t>(in.size());
        y /= static_cast<number_t>(in.size());

        return Point(x,y);
    }

    number_t stdDev(const PointContainer &in) {
        if(in.empty())
            return 0.0;
        DelaunayL2 T(in.begin(),in.end());
        using std::vector, std::cout;
        Point midPoint = midpoint(in);
        vector<number_t> distances;
        distances.reserve(in.size());

        for( auto p : in )
        {
            number_t dist = getDistance( p, midPoint );
            //cout<<dist<<"\n";
            distances.push_back(dist);
        }
        number_t meanDistanceFromMidpoint = std::accumulate( distances.begin(),
                                                       distances.end(),
                                                       0.0 ) / in.size();
        cout<<"mean:"<<meanDistanceFromMidpoint<<"\n";

        vector<number_t> squaredDifferences;
        squaredDifferences.reserve(in.size());

        for( auto dist : distances )
        {
            number_t diff = dist - meanDistanceFromMidpoint,
                    squaredDiff = pow(diff,2);
            squaredDifferences.push_back( squaredDiff );
        }

        number_t variance = std::accumulate( squaredDifferences.begin(),
                                       squaredDifferences.end(),
                                       0.0 ) / in.size(),
                standardDeviation = CGAL::sqrt( variance );

        cout<<"variance:"<<variance<<"\n";
        //cout<<"std-dev:"<<standardDeviation<<"\n";

        return standardDeviation;
    }

} //namespace spanner

#endif //LIBSPANNER_MEASURE_H
