/*
 * Wrappers and types to aid the use of CGAL's Linf Delaunay triangulation.
 *
 * Linf triangulation construction is slower than L2 but still much faster than TD.
 *
 * https://doc.cgal.org/latest/Segment_Delaunay_graph_Linf_2/index.html
 */

#ifndef LIBSPANNER_DELAUNAYLINF_H
#define LIBSPANNER_DELAUNAYLINF_H

#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

#include "../points/ordering.h"
#include "../types.h"

namespace spanner {


    typedef index_t info_t; // the type to store in VertexHandle->info()

    /*
     * CGAL Objects
     *
     * CGAL::Segment_Delaunay_graph_storage_traits_with_info_2
     * requires these functors, in their current form, although
     * our use of SDG does not result in these functors actually
     * being used by the SDG.
     *
     * They may cause "unused" warnings, but these are safe to
     * ignore.
     */

    template<typename T>
    struct InfoConvert {
        typedef T Info;
        typedef const Info &result_type;

        inline const Info &operator()(const Info &info0, [[maybe_unused]]bool) const {
            return info0; // just return the info of the supporting segment
        }

        inline const Info &operator()(const Info &info0, [[maybe_unused]]const Info &, [[maybe_unused]]bool) const {
            return info0; // just return the info of the supporting segment
        }
    };

    template<typename T>
    struct InfoMerge {
        typedef T Info;
        typedef Info result_type;

        inline Info operator()(const Info &info0, [[maybe_unused]]const Info &info1) const {
            return info0; // just return the info of the supporting segment
        }
    };

    typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<K, CGAL::Field_with_sqrt_tag> Gt;
    typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
            info_t,
        InfoConvert<info_t>,
        InfoMerge<info_t>> StorageTraits;

    typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt, StorageTraits> DelaunayLinf; // the main type


    void DelaunayLinfSpanner(const bdps::input_t &in, bdps::output_t &out) {

        const index_t n = in.size();
        if (n > SIZE_T_MAX - 1 || n <= 1) return;

        // construct Linf Delaunay triangulation
        bdps::input_t P(in);
        std::vector<size_t> index;
        spatialSort<K>(P, index);
        DelaunayLinf DT;
        DelaunayLinf::Site_2 site;
        index_t id = 0;

        // store the vertex handles
        std::vector<DelaunayLinf::Vertex_handle> handles(n);

        //FaceHandle hint;
        for (size_t entry : index) {
            Point p = P[entry];
            site = DelaunayLinf::Site_2::construct_site_2(p);
            auto vh = DT.insert(site, entry);
            //hint = vh->face();
            //vh->storage_site().info() = entry;
            handles[entry] = vh;
        }
        for (auto e = DT.finite_edges_begin(); e != DT.finite_edges_end(); ++e) {
            out.emplace_back(e->first->vertex((e->second + 1) % 3)->storage_site().info(),
                             e->first->vertex((e->second + 2) % 3)->storage_site().info());
        }
    }

    } // spanners

#endif //SPANNER_DELAUNAYLINF_H
