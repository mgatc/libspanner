//
// Created by matt on 2/2/22.
//

#ifndef LIBSPANNER_UTILITIES_H
#define LIBSPANNER_UTILITIES_H

#include <boost/functional/hash.hpp> // hashing pairs
#include <CGAL/point_generators_2.h>
#include <cstring>
#include <fstream>
#include <utility>

#include "constants.h"
#include "types.h"

namespace spanner {

    template<class T>
    bool contains(const T &V, const typename T::key_type &v) {
        return V.find(v) != V.end();
    }

    template<typename first_t, typename second_t>
    std::pair<second_t, first_t> reverse_pair(const std::pair<first_t, second_t> &in) {
        return std::make_pair(in.second, in.first);
    }

    template< typename T1, typename T2, typename F >
    void forBoth( const std::pair<T1,T2>& p, F func ) {
        func( p.first, p.second );
        func( p.second, p.first );
    }


/* If V contains v, remove v.
 * If V does not contain v, add it.
 * Return new value.
 */
//template< class T >
//bool toggle( T& V, const typename T::index_t& v ) {
//    bool found = contains( V, v );
//    if( found ) V.erase( v );
//    else V.insert( v );
//    return !found;
//}

    template<typename T>
    inline std::pair<T, T> makeNormalizedPair(const T &i, const T &j) {
        return std::make_pair(
                CGAL::min(i, j),
                CGAL::max(i, j)
        );
    }
    template< class T >
    inline std::pair<T, T> normalize_pair( const std::pair<T,T>& toNormalize ) {
        return makeNormalizedPair( toNormalize.first, toNormalize.second );
    }


    // FROM https://www.bfilipek.com/2018/06/variant.html#overload
    template<class... Ts> struct overload : Ts... { using Ts::operator()...; };
    template<class... Ts> overload(Ts...) -> overload<Ts...>;
    // END

    std::ostream& operator<<(std::ostream &os, mixed_t value) {
        visit( overload{
                [&](index_t& val ){ os<<val; },
                [&](number_t& val ){ os<<val; }
        }, value );
        return os;
    }
    std::string to_string(mixed_t value) {
        std::ostringstream oss;
        oss<<value;
        return oss.str();
    }


    // see https://stackoverflow.com/questions/5891610/how-to-remove-certain-characters-from-a-string-in-c
    std::string removeCharsFromString( std::string str, const char* charsToRemove ) {
        for ( unsigned int i = 0; i < std::strlen(charsToRemove); ++i ) {
            str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
        }
        return str;
    }
    std::string removeSpaces(std::string str) {
        const char space = ' ';
        return removeCharsFromString(str,&space);
    }
    std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
    {
        // https://stackoverflow.com/a/1120224
        std::vector<std::string>   result;
        std::string                line;
        std::getline(str,line);

        std::stringstream          lineStream(line);
        std::string                cell;

        while(std::getline(lineStream,cell, ','))
        {
            result.push_back(cell);
        }
        // This checks for a trailing comma with no data after it.
//        if (!lineStream && cell.empty())
//        {
//            // If there was a trailing comma then add an empty element.
//            result.push_back("");
//        }
        return result;
    }


    template<typename P>
    struct PointComparator {
        bool operator()(const P &p, const P &q) const noexcept {
            return CGAL::compare_xy(p,q) == CGAL::SMALLER;
        }
    };

    template< typename Point >
    struct PointHasher {
        std::size_t operator()(const Point& p) const noexcept {
            std::size_t seed = 31;
            boost::hash_combine( seed, p.x() );
            boost::hash_combine( seed, p.y() );
            return seed;
        }
    };

    template<typename T>
    struct MinHeapCompare {
        bool operator()( const T &n1, const T &n2 ) const {
            return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
        }
    };

    template<typename T>
    struct MaxHeapCompare {
        bool operator()( const T &n1, const T &n2 ) const {
            return (n1.first < n2.first) || ((n1.first == n2.first) && (n1.second < n2.second));
        }
    };

    struct IndexPairHash {
        size_t operator()(const index_tPair &item) const noexcept {
            size_t seed = 31;
            boost::hash_combine(seed, CGAL::min(item.first, item.second));
            boost::hash_combine(seed, CGAL::max(item.first, item.second));
            return seed;
        }
    };

    struct IndexPairComparator {
        bool operator()(const index_tPair &lhs, const index_tPair &rhs) const noexcept {
            return (lhs.first == rhs.first && lhs.second == rhs.second) ||
                   (lhs.first == rhs.second && lhs.second == rhs.first);
        }
    };

    struct PointConeHash {
        size_t operator()(const std::pair<index_t, cone_t> &PC) const noexcept {
            size_t seed = 31;
            boost::hash_combine(seed, PC.first);
            boost::hash_combine(seed, PC.second);
            return seed;
        }
    };

    struct PointConeComparator {
        bool operator()(const std::pair<index_t, cone_t> &lhs, const std::pair<index_t, cone_t> &rhs) const noexcept {
            return lhs.first == rhs.first && lhs.second == rhs.second;
        }
    };

    struct DegreeVertexComparator {
        bool operator()(const index_tPair &n1, const index_tPair &n2) const {
            return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
        }
    };



}

#endif //LIBSPANNER_UTILITIES_H
