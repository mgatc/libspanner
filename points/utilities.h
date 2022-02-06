#ifndef SPANNERS_POINTUTILITIES_H
#define SPANNERS_POINTUTILITIES_H

#include <fstream>
#include <iostream>
#include <vector>
//#include <cmath>
//#include <random>
//#include <string>
//#include <unordered_set>
//#include <utility>

#include "../types.h"
#include "../utilities.h"

namespace spanner {

    template<class Container, class Random>
    void perturb(Container &P, number_t val, Random rand = Random() ) {
        CGAL::perturb_points_2(P.begin(),P.end(),val,val,rand);
    }

    template< class OutputIterator >
    void readPointsFromFile( OutputIterator out, const std::string& outputFileName, const size_t n = SIZE_T_MAX) {
        std::ifstream in(outputFileName);
        if (in.is_open()) {
            number_t x,y;
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
    bool writePointsToFile(InputIterator begin, InputIterator end, std::string name="" ){
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
        return points.empty();
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

} // spanner

#endif //SPANNERS_POINTUTILITIES_H
