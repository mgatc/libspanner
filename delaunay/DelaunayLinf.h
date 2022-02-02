#ifndef GEOMETRIC_SPANNERS_DELAUNAYLINF_H
#define GEOMETRIC_SPANNERS_DELAUNAYLINF_H

#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

#include "types.h"

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
     * They will cause "unused" warnings, but these are safe to
     * ignore.
     */

    template<typename T>
    struct InfoConvert {
        typedef T Info;
        typedef const Info &result_type;

        inline const Info &operator()(const Info &info0, bool) const {
            return info0; // just return the info of the supporting segment
        }

        inline const Info &operator()(const Info &info0, const Info &, bool) const {
            return info0; // just return the info of the supporting segment
        }
    };

    template<typename T>
    struct InfoMerge {
        typedef T Info;
        typedef Info result_type;

        inline Info operator()(const Info &info0, const Info &info1) const {
            return info0; // just return the info of the supporting segment
        }
    };

    typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<K, CGAL::Field_with_sqrt_tag> Gt;
    typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
            info_t,
        InfoConvert<info_t>,
        InfoMerge<info_t>> StorageTraits;

    typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt, StorageTraits> DelaunayLinf; // the main type

} // spanners

#endif //GEOMETRIC_SPANNERS_DELAUNAYLINF_H
