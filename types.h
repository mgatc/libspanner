#ifndef LIBSPANNER_TYPES_H
#define LIBSPANNER_TYPES_H

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <boost/functional/hash.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace spanner {

    // The geometry kernel aka how to represent the math
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;

    typedef Epick K;

    typedef K::Point_2 Point;
    typedef std::vector<Point> PointContainer;
    typedef K::Vector_2 Vector_2;

    typedef K::FT number_t;
    typedef size_t index_t;
    typedef size_t cone_t;

    typedef std::variant<index_t,number_t> mixed_t;

    typedef std::pair<index_t, index_t> index_tPair;
    typedef index_tPair Edge;
    typedef boost::hash<index_tPair> index_tPairHash;
    typedef std::unordered_set<index_tPair, index_tPairHash> index_tPairSet;
    typedef std::unordered_map<index_tPair, bool, index_tPairHash> index_tPairMap;

    typedef std::vector<std::unordered_set<index_t>> AdjacencyListDense;
    typedef std::unordered_map<index_t,std::unordered_set<index_t>> AdjacencyListSparse;

}

#endif //LIBLIBSPANNER_TYPES_H
