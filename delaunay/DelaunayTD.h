/*
 * Wrappers and types to aid the use of CGAL's Theta graph construction to create a TD Delaunay triangulation.
 *
 * TD Delaunay construction is quite slow.
 *
 * https://doc.cgal.org/latest/Cone_spanners_2/classCGAL_1_1Construct__theta__graph__2.html
 *
 * uses boost graph library (bgl) to represent the graph
 * https://www.boost.org/doc/libs/1_78_0/libs/graph/doc/
 */

#ifndef LIBSPANNER_DELAUNAYTD_H
#define LIBSPANNER_DELAUNAYTD_H

#include <CGAL/Compact_container.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/property_map.h>

#include <boost/functional/hash.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include "../constants.h"
#include "../geometry.h"
#include "../types.h"
#include "../utilities.h"

namespace spanner {

    namespace td {
//Cone angles.
        const number_t alpha = PI / 3;

//Slopes of the cone boundary lines.
        const number_t orthBisectorSlopes[] = {0, -1 * COT30, COT30, 0, -1 * COT30, COT30};

//Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an getAngle of (PI/3).
        cone_t getSingleCone(const index_t p, const index_t q, const std::vector <Point> &H) {
            if (CGAL::compare_y(H[p], H[q]) == CGAL::EQUAL) {
                return 1 + 3 * int(CGAL::compare_x(H[p], H[q]) == CGAL::LARGER);
            }
            const Point refPoint(H.at(p).x() - TAN30, H[p].y() + 1);
            const number_t theta = getAngle(refPoint, H[p], H.at(q));
            const cone_t cone = (theta / alpha);

            return cone;
        }

//Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
        size_t getCone(const size_t p, const size_t q, const std::vector <Point> &H) {
            return p < q ?
                   getSingleCone(p, q, H)
                         : (getSingleCone(q, p, H) + 3) % 6;
        }

        template<class Gt>
        class HalfThetaTriangulation {
        public:

            // select the kernel type
            typedef Gt Geom_traits;
            typedef typename Geom_traits::Point_2 Point_2;
            typedef typename Geom_traits::Direction_2 Direction_2;
            // typedef Vertex_base<Geom_traits>                        VertexHandle;
            typedef typename std::vector<Point_2>::iterator Point_iterator;
            //typedef typename std::vector<VertexHandle>::iterator   Vertex_iterator;
            typedef boost::adjacency_list<boost::vecS,
                    boost::vecS,
                    boost::bidirectionalS,
                    Point_2,
                    size_t> Graph;
            typedef typename boost::graph_traits<Graph>::edge_iterator Edge_iterator;
            typedef typename boost::graph_traits<Graph>::in_edge_iterator In_edge_iterator;
            typedef typename boost::graph_traits<Graph>::out_edge_iterator Out_edge_iterator;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_descriptor;
            typedef typename boost::graph_traits<Graph>::vertex_iterator Vertex_iterator;
            typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;

            // define the Graph

            Graph _G;
            CGAL::Construct_theta_graph_2<Geom_traits, Graph> halfTheta;

            template<typename InputIt>
            HalfThetaTriangulation(InputIt first, InputIt last, const Geom_traits &gt = Geom_traits()) :
                    halfTheta(6, Direction_2(1, 0), CGAL::ODD_CONES),
                    _gt(gt) {
                insert(first, last);
            }

//    HalfThetaTriangulation(const Geom_traits& gt = Geom_traits()) :
//                halfTheta(6, Direction_2(1,0), CGAL::ODD_CONES),
//                _gt(gt) {}

            template<class InputIt>
            std::ptrdiff_t insert(InputIt first, InputIt last) {
                auto n = number_of_vertices();

                _P.insert(_P.end(), first, last);
                size_t id = 0;
                for (auto pit = _P.begin(); pit != _P.end(); pit++) {
                    info.emplace(pit, id++);
                }
                halfTheta(_P.begin(), _P.end(), _G);          // Construct the half-theta graph


                return number_of_vertices() - n;
            }

//    const Geom_traits& geom_traits() const { return _gt;}

            std::ptrdiff_t number_of_vertices() {
                return boost::num_vertices(_G);
            }

            Vertex_iterator finite_vertices_begin() const {
                return boost::vertices(_G).first;
            }

            // finite_vertices_end
            Vertex_iterator finite_vertices_end() const {
                return boost::vertices(_G).second;
            }

            // finite_points_begin
            Point_iterator points_begin() {
                return _P.begin();
            }

            // finite_points_end
            Point_iterator points_end() {
                return _P.end();
            }

            // finite_edges_begin
            Edge_iterator edges_begin() const {
                return boost::edges(_G).first;
            }

            // finite_edges_end
            Edge_iterator edges_end() const {
                return boost::edges(_G).second;
            }

            Point_2 &point(size_t index) {
                return _G[index];
            }

            // get negative cone edges
            In_edge_iterator negative_cone_edges_begin(const VertexDescriptor &u) const {
                return in_edges(u, _G).first;
            }

            In_edge_iterator negative_cone_edges_end(const VertexDescriptor &u) const {
                return in_edges(u, _G).second;
            }

            // get positive cone edges
            Out_edge_iterator positive_cone_edges_begin(const VertexDescriptor &u) const {
                return out_edges(u, _G).first;
            }

            Out_edge_iterator positive_cone_edges_end(const VertexDescriptor &u) const {
                return out_edges(u, _G).second;
            }

            //Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
            inline size_t getCone(const VertexDescriptor &p, const VertexDescriptor &q) const {
                return td::getCone(p, q, _P);
            }

            inline bool edgeExists(const std::pair <VertexDescriptor, VertexDescriptor> &e) const {
                return boost::edge(e.first, e.second, _G).second;
            }

            inline std::pair<std::pair<VertexDescriptor, VertexDescriptor>, bool>
            eitherEdge(const VertexDescriptor &u, const VertexDescriptor &v) const {
                auto orientedEdge = std::make_pair(u, v),
                        reversedEdge = reverse_pair(orientedEdge);
                bool reversedEdgeExists = edgeExists(reversedEdge);
                return make_pair(
                        reversedEdgeExists ? reversedEdge : orientedEdge,
                        reversedEdgeExists || edgeExists(orientedEdge)
                );
            }

            template<class EdgeList>
            void fanOfCone(const VertexDescriptor &v, const size_t cone, EdgeList &fan) const {
                // get in_edges of v
                for (auto neit = negative_cone_edges_begin(v);
                     neit != negative_cone_edges_end(v); ++neit) {
                    auto e = *neit;
                    if (getCone(target(e), source(e)) == cone) {
                        fan.insert(fan.end(), e);
                    }
                }
                // sort edges in counter-clockwise order
                sort(fan.begin(), fan.end(),
                     [&](const auto &lhs, const auto &rhs) {
                         //assert( target(lhs) == target(rhs) );
                         const Point_2 refPoint(_P[target(lhs)].x() - TAN30,
                                                _P[target(lhs)].y() + 1);
                         return getAngle(refPoint,
                                         _P[target(lhs)],
                                         _P[source(lhs)])
                                > getAngle(refPoint,
                                           _P[target(rhs)],
                                           _P[source(rhs)]);
                     }
                );
//        cout<< "Fan "<<cone<< " of "<<v<<": ";
//        for( auto e : fan )
//        {
//            cout<<source(e)<<"-"<<target(e)<<" ";
//        }
//        cout<<"\n";
            }

            VertexDescriptor source(const Edge_descriptor &e) const {
                return boost::source(e, _G);
            }

            VertexDescriptor target(const Edge_descriptor &e) const {
                return boost::target(e, _G);
            }

            VertexDescriptor parent(const VertexDescriptor child, const size_t i) const {
                for (auto it = positive_cone_edges_begin(child);
                     it != positive_cone_edges_end(child); ++it) {
                    VertexDescriptor v = target(*it);
                    if (getCone(child, v) == i)
                        return v;
                }
                return SIZE_T_MAX;
            }


        protected:
            Geom_traits _gt;
            //VertexHandle _infinite_vertex;
            std::vector<Point_2> _P;
            std::vector<VertexDescriptor> _V;
            std::map<Point_iterator, size_t> info;

        private:
            const double alpha = PI / 3;
        };
    }

    typedef td::HalfThetaTriangulation<K> DelaunayTD; // the main type


    void DelaunayTDSpanner(const bdps::input_t &in, bdps::output_t &out) {

        const index_t n = in.size();
        if (n > SIZE_T_MAX - 1 || n <= 1) return;

        std::vector<Point> P(in);

        DelaunayTD D(P.begin(), P.end());

        std::set<Edge> E;

        for(auto e=D.edges_begin(); e!=D.edges_end(); ++e) {
            E.emplace(std::min(D.source(*e),D.target(*e)), std::max(D.source(*e), D.target(*e)));
        }
        std::copy(E.begin(), E.end(), std::back_inserter(out));
    }
}

#endif //TD_DELAUNAY_H
