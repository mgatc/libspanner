//
// Created by matt on 7/15/21.
//

#ifndef LIBSPANNER_CLEANER_H
#define LIBSPANNER_CLEANER_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/property_map.h>

#include "../utilities.h"

namespace spanner {

template <class Geom_traits>
class PointSetCleaner {
public:

    typedef std::pair<std::string,std::string> PointString;
    std::vector<PointString> P;

    //typedef Geom_traits Gt;
    typedef typename Geom_traits::Point_2 Point_2;
    //typedef typename Geom_traits::FT FT;

    bool addFromFile(const std::string& filename) {
        using std::cout;
        using CGAL::IO::read_XYZ;

        const size_t n_prev = this->size();

        std::ifstream in(filename);
        if (!in.is_open())
            assert(!"ERROR WRITING FILE");

        CGAL::IO::set_ascii_mode(in);

        std::vector<typename Geom_traits::Point_3> temp3d;


//        if(!read_XYZ(in, back_inserter(temp3d))) {
//            assert(!"Error opening file!\n");
//        }
//
//        std::transform(temp3d.begin(), temp3d.end(), std::back_inserter(*this),
//            [] (const auto& p) {
//                return typename Geom_traits::Point_2{ p.x(), p.y() };
//            });

        Point_2 p;
        while (in >> p) {
            this->push_back(p);
        }
        in.close();

        cout << "Added " << (this->size() - n_prev) << " points from file\n";

        return true;
    }

//    bool writeToFile(const std::string& filename) {
//        if(P.empty()) {
//            std::cerr<< "Point set is empty\n";
//            return false;
//        }
//
//        std::ofstream out(filename);
//        CGAL::IO::set_ascii_mode(out);
//
//        bool result = CGAL::IO::write_XYZ(out,*this);
//        if(!result) std::cerr << "ERROR WRITING TO FILE "<<filename<<std::endl;
//        return result;
//    }

    size_t removeDuplicates() {
        const size_t n_prev = P.size();

        for(auto p : P) {
            std::cout<<p.first<<" "<<p.second<<"\n";
        }
        typedef typename Geom_traits::Point_2 Point;
        std::map<Point, PointString, spanner::PointComparator<Point>> Temp;
        std::transform(P.begin(),P.end(),std::inserter(Temp, Temp.end()),[](const auto& p) {
            typename Geom_traits::FT x, y;
            std::stringstream str;
            str << p.first;
            str >> x;
            str << p.second;
            str >> y;

            return std::make_pair(Point(x,y), p);
        });
        P.clear();
        std::transform(Temp.begin(),Temp.end(),std::back_inserter(P),[](const auto& p) {
            return p.second;
        });

        return n_prev - P.size();
    }

    bool writeToFile(const std::string& filename) {
        std::ofstream out;
//        string fName;
//        fName = "data-" + to_string(n) + "_" + to_string(size) + "x" + to_string(size) + ".txt";
        out.open(filename, std::ios::trunc);

        if (out.is_open()) {
            for (auto it=P.begin(); it!=P.end(); ++it) {
                std::string p(it->first + " " + it->second);
                out << p << "\n";
            }
            out.close();
            return true;
        } else {
            std::cout << "ERROR WRITING TO FILE " << filename <<std::endl;
            return false;
        }
    }

    PointSetCleaner(const std::string& filename,
                       const size_t xCol,
                       const size_t yCol,
                       const char delimiter = ',',
                       const std::string& charsToRemove = "" ) {
        using std::cout, std::endl;
        using std::string;

        typedef typename Geom_traits::FT number_t;

        std::ifstream in(filename);
        assert(in.is_open() && "Failed to open input file!");

        size_t discarded = 0;

        string line;
        string cell;
        auto chars = charsToRemove.c_str();

        while (getline(in, line)) {
            std::stringstream linestream(line);
            std::vector<string> result;
            while (getline(linestream, cell, delimiter)) {
                cell = spanner::removeCharsFromString( cell, chars );
                if (result.size() == xCol)
                    cout << "x) " << cell << "\n";
                if (result.size() == yCol)
                    cout << "y) " << cell << "\n";
                result.push_back(cell);
            }
            cout << "\n";

            bool success = result.size() > std::max(xCol, yCol);
            if (success) {
                try {
                    std::stod(result.at(xCol)); //test that it passes stod
                    std::stod(result.at(yCol));
                    P.emplace_back(result.at(xCol), result.at(yCol)); // store original string from file
                }
                catch (std::invalid_argument &ia) {
                    success = false;
                }
            }
            if (!success) {

                cout << "DISCARDED ROW: " << line << "\n";
                //cout<< "\""<<result.at(xCol)<<"\", \""<< result.at(yCol)<<"\""<<std::endl;
                //cout<< "line "<<(P.size()-1)<<endl<<endl;
                cout<< "Continue? Y/n"<<endl;
                char response = 'x';
                while(response != 'y' && response !='n' ) {
                    std::cin>>response;
                    response = std::tolower(response);
                }
                if(response == 'n')
                    assert(!"discarded row");

                cout << "DISCARDED ROW: " << line << "\n";

                ++discarded;
            }
        }
        in.close();

        string tempFileOut (filename + ".temp"),
                finalFileOut(filename + ".txt");

        cout << "Found "<<P.size()<<" points\n"<<std::endl;

        auto dupes = removeDuplicates();
        cout<< "Removed "<<dupes<<" duplicates\n"<<std::endl;
        cout<< "Final count: "<<P.size()<<endl;
        writeToFile(finalFileOut);
    }
};

} // namespace spanner

#endif //LIBSPANNER_CLEANER_H
