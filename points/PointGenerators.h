#ifndef SPANNERS_POINTGENERATORS_H
#define SPANNERS_POINTGENERATORS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>

#include <boost/functional/hash.hpp>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include "tools/Utilities.h"

namespace bdps_experiment {

    enum DistributionType {
        DistributionTypeFirst = 0,
        Synthetic = DistributionTypeFirst,
        Real,
        DistributionTypeLast
    };

    enum SyntheticDistribution {
        SyntheticDistributionFirst=0,
//        UniformInsideSquare,
//        UniformInsideDisc,
////        UniformOnSquare,
////        UniformOnCircle,
//        NormalInsideSquare,
        NormalClustersInsideSquare = SyntheticDistributionFirst,
//        ContiguousGrid,
//        UniformRandomGrid,
//        UniformInsideAnnulus,
//        Galaxy,
//        ConvexHullInDisc,// = SyntheticDistributionFirst,
        SyntheticDistributionLast
    };

    const std::vector<std::string> SYNTHETIC_DISTRIBUTION_NAMES = {
//            "Uniform Inside Square",
//            "Uniform Inside Disc",
////            "Uniform On Square",
////            "Uniform On Circle",
//            "Normal Inside Square",
            "Normal Inside Square with Clusters",
//            "Contiguous Grid",
//            "Uniform Random Grid",
//            "Uniform Inside Annulus",
//            "Galaxy",
//            "Convex Hull In Disc"
    };
    std::vector<std::string> REAL_POINTSET_NAMES;

    class PointGenerator_2 {

    public:
        PointGenerator_2() : m_randCgal(std::rand()) {}

        void loadFromFile(const std::string filename, std::vector<Point> &P) {
            std::ifstream in(filename);

            if (!in.is_open())
                cout<<"Error opening file!\n";

            number_t x,y;
            while ( in >> x >> y ) {
                P.emplace_back(x,y);
            }
            in.close();

            perturb(P, m_perturbationValue);
        }

        void insideSquare(const index_t n, const double sizeOfSquare, std::vector<Point> &P) {
            typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            std::unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }
        void onSquare(const index_t n, const double sizeOfSquare, std::vector<Point> &P) {
            typedef CGAL::Random_points_on_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            std::unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }
        void onCircle(const index_t n, const double sizeOfSquare, std::vector<Point> &P) {
            typedef CGAL::Random_points_on_circle_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            std::unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }

        void onSpokes(const index_t n, const unsigned numSpokes, std::vector<Point> &P) {
            //srand(seed());
            double spokeAngle = 2*PI / numSpokes;

            std::unordered_set<Point> P_unique;

            while( P_unique.size() < n ) {
                double distance = randFloat();
                double angle = randFloat() * 2 * PI;
                angle = ((unsigned) (angle / spokeAngle) )*spokeAngle;
                P_unique.emplace(cos(angle) * distance,
                                 sin(angle) * distance);
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }
        void inGalaxy(const index_t n, const unsigned numSpokes, std::vector<Point> &P) {
            // see https://itinerantgames.tumblr.com/post/78592276402/a-2d-procedural-galaxy-with-c
            //srand(seed());
            const double spokeAngle = 2*PI / numSpokes,
                    armOffsetMax = 0.5,
                    rotationFactor = 5,
                    perturbationValue = 0.02;

            std::unordered_set<Point> P_unique;

            while( P_unique.size() < n ) {
                //for(index_t i=0; i<n; ++i) {
                double distance = randFloat();
                distance = pow(distance,2);

                double angle = randFloat() * 2 * PI;
                double armOffset = randFloat() * armOffsetMax;
                armOffset -= armOffsetMax / 2;
                armOffset *= (1/distance);

                double squaredArmOffset = pow(armOffset,2);
                squaredArmOffset *= -1 * int(armOffset < 0);
                armOffset = squaredArmOffset;

                double rotation = distance * rotationFactor;

                angle = ((unsigned) (angle / spokeAngle) )*spokeAngle;
                angle += armOffset + rotation;

                P_unique.emplace(cos(angle) * distance,
                                 sin(angle) * distance);
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, perturbationValue);
        }

        void insideDisc(const index_t n, const double radiusOfDisk, std::vector<Point> &P) {
            typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(radiusOfDisk, m_randCgal);

            std::unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void onConvexHullInDisc(const index_t n,
                                const double radius,
                                std::vector<Point> &P) {
            typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            std::unordered_set<Point> P_unique;
            size_t remaining;

            while((remaining = n - P_unique.size()) > 0)
                CGAL::random_convex_set_2(remaining,inserter(P_unique, P_unique.end()),Point_generator(radius, m_randCgal));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }


        void insideSquareNormal(const index_t pointsInACuster,
                                const index_t numberOfClusters,
                                std::vector<Point> &P,
                                const number_t xStdDev = 2.0,
                                const number_t yStdDev = 2.0) {
            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::default_random_engine generatorX(rngX()), generatorY(rngY());
            std::normal_distribution<double> distributionX(0.0, xStdDev), distributionY(2.0, yStdDev);

            std::mt19937 rngShift(std::random_device{}());

            std::uniform_int_distribution shiftDistribution(0, INT32_MAX);

            index_t shiftX, shiftY;
            std::unordered_set<std::pair<index_t, index_t>, boost::hash<std::pair<index_t, index_t>>> S;

            std::unordered_set<Point> P_unique;

            for (index_t c = 0; c < numberOfClusters; c++) {
                if (c != 0) {
                    shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                    shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);

                    while (!(S.find(make_pair(shiftX, shiftY)) == S.end())) {
                        shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                        shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);
                    }
                } else
                    shiftX = shiftY = 0;

                S.insert(make_pair(shiftX, shiftY));

                for (index_t i = 0; i < pointsInACuster; i++) {
                    double x = distributionX(generatorX) + shiftX;
                    double y = distributionY(generatorY) + shiftY;
                    P_unique.emplace(x, y);
                }
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void contiguousOnGrid(const index_t n, std::vector<Point> &P) {
            CGAL::points_on_square_grid_2(ceil(std::sqrt(n)), n, std::back_inserter(P), CGAL::Creator_uniform_2<number_t, Point>());
            perturb(P, m_perturbationValue);
        }

        void randomOnGrid(const index_t n, std::vector<Point> &P) {
            std::unordered_set<pair<int, int>, boost::hash<pair<int, int>>> S;
            std::unordered_set<Point> P_unique;

            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::uniform_int_distribution xDistribution(0, (int) ceil(0.7 * n)), yDistribution(0, (int) ceil(0.7 * n));

            index_t count = 0;

            while (count < n) {
                int x = xDistribution(rngX), y = yDistribution(rngY);

                if (S.find(make_pair(x, y)) == S.end()) {
                    P_unique.emplace(x, y);
                    S.insert(make_pair(x, y));
                    count++;
                }
            }

            std::copy(P_unique.begin(),P_unique.end(),std::back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void insideAnnulus(const index_t n, const double r2, const double r1, std::vector<Point> &P) {
            assert(r2 > r1);
            std::unordered_set<Point> P_unique;

            std::default_random_engine generator(seed());
            std::uniform_real_distribution<double> distributionR(r1, r2), distributionT(0, 1);

            for (index_t i = 0; i < n; i++) {
                double t = 2 * M_PI * distributionT(generator);
                double r = distributionR(generator);
                P_unique.emplace(r * cos(t), r * sin(t));
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

    private:
        CGAL::Random m_randCgal;
        inline static number_t m_perturbationValue = 0.0001;

        size_t seed() {
            return std::rand();
        }
        double randFloat() {
            return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        }
        template<class Container>
        void perturb(Container &P, number_t val) {
            CGAL::perturb_points_2(P.begin(), P.end(), val, val,m_randCgal);
        }
    };



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
        return points.empty();
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

        index_t remaining;
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


} // bdps_experiment

#endif //SPANNERS_POINTGENERATORS_H
