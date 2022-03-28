#ifndef LIBSPANNER_GREEDYSPANNERS_H
#define LIBSPANNER_GREEDYSPANNERS_H

#include <string>
#include <vector>

#include "greedy/DN97.h"

namespace spanner {

    enum GreedySpannerAlgorithm {
        GreedySpannerAlgorithmFirst = 0,
        Dn97 = AlgorithmFirst,
        GreedySpannerAlgorithmLast
    };

    namespace greedy {

        const std::string ALGORITHM_SYMBOL = "Algorithm";
        const std::vector<std::string> ALGORITHM_NAMES = {
                "DN97",
        };
    } // greedy

} // spanner

#endif //LIBSPANNER_GREEDYSPANNERS_H
