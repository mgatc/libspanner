//
// Created by matt on 2/2/22.
//

#ifndef LIBSPANNER_UTILITIES_H
#define LIBSPANNER_UTILITIES_H

#include <boost/functional/hash.hpp> // hashing pairs
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
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

    template< typename Point >
    struct PointHasher {
        std::size_t operator()(const Point& p) const noexcept {
            std::size_t seed = 31;
            boost::hash_combine( seed, p.x() );
            boost::hash_combine( seed, p.y() );
            return seed;
        }
    };

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
        return make_pair(
                CGAL::min(i, j),
                CGAL::max(i, j)
        );
    }
    template< class T >
    inline std::pair<T, T> normalize_pair( const std::pair<T,T>& toNormalize ) {
        return makeNormalizedPair( toNormalize.first, toNormalize.second );
    }

    template<typename Point_2>
    inline number_t getDistance(const Point_2 &p, const Point_2 &q) {
        return p == q ? 0 : CGAL::sqrt(CGAL::squared_distance(p, q));
    }

    template< class Points >
    inline number_t getEdgeLength(const Edge& e, const Points& P ) {
        return getDistance(P[e.first], P[e.second]);
    }

    template<class Point_2>
    inline number_t getAngle(const Point_2 &p, const Point_2 &q, const Point_2 &r) {
        auto pq = q - p,
                rq = q - r;

        number_t result = atan2(pq.y(), pq.x()) - atan2(rq.y(), rq.x());

        // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
        // Our zero is also "up," but we only want positive values between 0 and 2*PI:

        result = fmod(result + 2 * PI, 2 * PI);

//        cout<<p<<" - "<<q<<" - "<<r<<endl;
//        cout<<"angle="<<result<<endl;
        return result;
    }

    template<class Points>
    inline number_t getAngle(const index_t p, const index_t q, const index_t r, const Points &P) {
        return getAngle(P[p], P[q], P[r]);
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
        oss<<std::move(value);
        return oss.str();
    }


    std::string removeSpaces(std::string str) {
        str.erase(std::remove_if(str.begin(), str.end(),
                                 [](auto x) { return std::isspace(x); }), str.end());
        return str;
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


    template< class OutputIterator >
    void readPointsFromFile( OutputIterator out, const std::string& outputFileName, const size_t n=SIZE_T_MAX ) {
        std::ifstream in(outputFileName);
        if (in.is_open()) {
            double x,y;
            size_t i = 0;
            while ( i<n && in >> x >> y ) {
                *out = Point(x,y);
                ++out;
                ++i;
            }
            in.close();
        }
    }

    template<class InputIterator>
    bool writePointsToFile(InputIterator begin, InputIterator end, std::string name="") {
        std::vector<Point> points(begin,end);
        std::ofstream out;
        if(name.empty())
            name = "data-" + to_string(points.size()) + ".xy";
        out.open( name, std::ios::trunc );

        if(!out.is_open())
            return false;

        for( Point p : points )
            out << p << std::endl;

        out.close();
        return points.size() > 0;
    }

    template< class OutputIterator >
    void generateRandomPoints( index_t n, number_t size, OutputIterator pointsOut ) {
        typedef CGAL::Creator_uniform_2<number_t,Point> Creator;

        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( size );
        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   size );
        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );


        auto g1s = CGAL::Random_points_in_square_2<Point,Creator>( size/4 );
        auto g2s = CGAL::Random_points_in_disc_2<Point,Creator>(   size/4 );
        auto g3s = CGAL::Random_points_on_square_2<Point,Creator>( size/4 );
        auto g4s = CGAL::Random_points_on_circle_2<Point,Creator>( size/4 );

        std::set<Point> points;

//        std::copy_n( g2, n/9, inserter(points) );
//        std::copy_n( g3, n/9, inserter(points) );
//        std::copy_n( g4, n/9, inserter(points) );
//
//        std::copy_n( g1s, n/9, inserter(points) );
//        std::copy_n( g2s, n*2/9, inserter(points) );
//        std::copy_n( g3s, n/9, inserter(points) );
//        std::copy_n( g4s, n/18, inserter(points) );

        int remaining;
        while( (remaining = n - points.size()) > 0 ) {
            std::copy_n( g1, remaining, inserter(points) );
        }


        //points.emplace(0,0);

        // copy points to output iterator
        for( Point p : points )
            *(pointsOut++) = p;

        writePointsToFile(points.begin(),points.end());

    }


//    template< class OutputIterator >
//    string generatePointsNovel( OutputIterator pointsOut, size_t rows = 10, size_t cols = 10 ) {
//
//        vector<Point> points;
//
//        const double skew = 0.01;
//        for( size_t i=0; i<rows; ++i ) {
//            bool rowIsOdd = i%2;
//            for( size_t j=rowIsOdd; j<cols; j+=1+(rowIsOdd) ) {
//                bool colIsOdd = j%2;
//                double y = static_cast<double>(i);
//                y += (rowIsOdd || colIsOdd) ?
//                     0 : (skew * ( (i+j)%4 == 0 ? -1 : 1 ) );
////            if( rowIsEven && j%2 == 0 ) {
////                if( (i+j)%4 == 0 ) {
////                    y -= skew;
////                } else {
////                    y += skew;
////                }
////            }
//                Point p(j,y);
//                //cout<<p<<"\n";
//                points.push_back(p);
//            }
//        }
//
//        // copy points to output iterator
//        for( Point p : points )
//            *(pointsOut++) = p;
//
//        // copy points to file
//        ofstream out;
//        string fName;
//        fName = "data-NOVEL-" + to_string(points.size()) + "_" + to_string(rows) + "x" + to_string(cols) + ".txt";
//        out.open( fName, ios::trunc );
//        for( Point p : points )
//            out << p << endl;
//
//        out.close();
//
//        return fName;
//    }

    template< typename Point_2 >
    struct EuclideanDistanceToPoint {
        Point_2 goal;
        double operator()( Point_2 p ) {
            return getDistance(p, goal);
        }
    };

    template<typename T>
    struct MinHeapCompare {
        bool operator()( const T &n1, const T &n2 ) const {
            return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
        }
    };

//template<typename T>
//struct MaxHeapCompare {
//    bool operator()( const T &n1, const T &n2 ) const {
//        return (n1.first < n2.first) || ((n1.first == n2.first) && (n1.second < n2.second));
//    }
//};

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

    template<class K>
    void spatialSort(std::vector<typename K::Point_2> &P, std::vector<index_t> &index) {
        typedef CGAL::Spatial_sort_traits_adapter_2<K,
                typename CGAL::Pointer_property_map<typename K::Point_2>::type> SearchTraits;

        index.clear();
        index.reserve(P.size());

        std::copy(boost::counting_iterator<std::size_t>(0),
                  boost::counting_iterator<std::size_t>(P.size()),
                  std::back_inserter(index));

        CGAL::spatial_sort(index.begin(),
                           index.end(),
                           SearchTraits(CGAL::make_property_map(P)));
        //cout<<"done sorting"<<endl;
    }


    struct DegreeVertexComparator {
        bool operator()(const index_tPair &n1, const index_tPair &n2) const {
            return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
        }
    };

    template<class Graph, class OutputIterator>
    void reverseLowDegreeOrdering(const Graph &T, OutputIterator out) {

        typedef boost::heap::fibonacci_heap<index_tPair, boost::heap::compare<DegreeVertexComparator>> Heap;
        typedef Heap::handle_type HeapHandle;
        typedef typename Graph::Vertex_circulator VertexCirculator;

        const index_t n = T.number_of_vertices();

        Heap H;
        vector<HeapHandle> handleToHeap(n);
        //vector<size_t> piIndexedByV(n);
        vector<index_t> ordering(n);
        vector<unordered_set<index_t>> currentNeighbors(n);

        // Initialize the vector currentNeighbors with appropriate neighbors for every vertex
        for (auto it = T.finite_vertices_begin();
             it != T.finite_vertices_end(); ++it) {
            VertexCirculator N = T.incident_vertices(it),
                    done(N);
            do {
                if (!T.is_infinite(N))
                    currentNeighbors.at(it->info()).insert(N->info());
            } while (++N != done);

            size_t degree = currentNeighbors.at(it->info()).size();
            handleToHeap[it->info()] = H.emplace(degree, it->info());
        }

        // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
        index_t i = n - 1; // start at the last valid index

        while (!H.empty()) {
            index_tPair p = H.top();
            H.pop();
            // make sure our math is correct, e.g., degree from heap key matches neighbor container size
            //assert(p.first == currentNeighbors.at(p.second).size());
            //assert(0 <= p.first && p.first <= 5); // Lemma 1

            // Erase this vertex from incidence list of neighbors and update the neighbors' key in the heap
            for (index_t neighbor : currentNeighbors.at(p.second)) {
                currentNeighbors.at(neighbor).erase(p.second);
                HeapHandle h = handleToHeap.at(neighbor);
                index_tPair q = make_pair(currentNeighbors.at(neighbor).size(), neighbor);
                H.update(h, q);
                H.update(h);
            }
            currentNeighbors.at(p.second).clear();
            //piIndexedByV[p.second] = i;
            ordering[i] = p.second;
            --i;
        }
        std::copy(ordering.begin(), ordering.end(), out );

//        for (auto v : ordering) {
//            *out = v;
//            ++out;
//        }
    }
    template< class DelaunayTriangulation, class VertexHash, class VertexHandle>
    inline index_t incidentChords(const DelaunayTriangulation& DT, const VertexHash& onOuterFace, const VertexHandle& v_k ) {
        typedef typename DelaunayTriangulation::Vertex_circulator VertexCirculator;
        VertexCirculator N = DT.incident_vertices(v_k),
                done(N);
        size_t c = 0;
        do {
            if( contains(onOuterFace, N->handle() ) )
                ++c;
        } while( ++N != done );
        // v_k is guaranteed to be incident to two vertices in onOuterFace, its neighbors.
        // Anything >2 is a chord
        //assert( c >= 2 );
        return (c - 2);
    }


    /**
     *  Given a Delaunay Triangulation DT and an output list out, compute the canonical out of
     *  the underlying point set.
     */
    template< class DelaunayTriangulation, typename OutputIterator >
    void canonicalOrder( const DelaunayTriangulation& DT, OutputIterator out ) {
        //Timer t(",");
        typedef typename DelaunayTriangulation::Vertex_handle VertexHandle;
        typedef typename DelaunayTriangulation::Vertex_circulator VertexCirculator;
        typedef unordered_set<VertexHandle> VertexHash;

        VertexHash onOuterFace, complete;
        queue<VertexHandle> ready;
        index_t i = DT.number_of_vertices();

        vector<VertexHandle> ordering(i);

        VertexCirculator v_convexHull = DT.incident_vertices(DT.infinite_vertex() ), // create a circulator of the convex hull
        done(v_convexHull );

        do { // cycle through convex hull to set vertex info
            onOuterFace.insert(v_convexHull->handle() );
            ready.push(v_convexHull->handle() );
        } while(++v_convexHull != done );

        // Reserve v_1 and v_2 so we can guarantee they are on the convex hull
        for(index_t j=0; j < 2; ++j ) {
            ordering[j] = ready.front();
            ready.pop();
        }

        VertexHandle v_k = ready.front();

        while( !ready.empty() ) {
            v_k = ready.front();
            //std::cout<<v_k->point()<<" ";
            ready.pop();
//            _algoTV.addToEventQueue( v_k, 0 );
//            _algoTV.addToEventQueue( v_k, false );

            if(incidentChords(DT, onOuterFace, v_k ) > 0 ) {
                ready.push(v_k);
                //std::cout<<"requeued";
            } else {
                //std::cout<<"processed";
                ordering[--i] = v_k;
                onOuterFace.erase(v_k);
                // add all neighbors not Complete or on outer face to ready list
                // add all neighbors not Complete to on outer face
                VertexCirculator N = DT.incident_vertices(v_k);
                done = N;
                do {
                    if( !contains( complete, N->handle() ) && !DT.is_infinite( N->handle() ) ) {
                        if( !contains(onOuterFace, N->handle() ) )
                            ready.push( N->handle() );
                        onOuterFace.insert(N->handle() );
                    }
                } while( ++N != done );
                complete.insert(v_k);
            }
            //std::cout<<"\n";
        }
        std::copy(ordering.begin(), ordering.end(), out );
    }
}

#endif //LIBSPANNER_UTILITIES_H
