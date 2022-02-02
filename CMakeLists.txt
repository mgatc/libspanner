cmake_minimum_required(VERSION 3.19)
project(libspanner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

## UNCOMMENT FOR RELEASE
#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")
#set(CMAKE_BUILD_TYPE Release)

include_directories(/usr/include)
include_directories(include)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include)

include_directories(.)
include_directories(bdps)

add_executable(libspanner

        bdps/BCC2012.h
        bdps/BGHP2010.h
        bdps/BGS2005.h
        bdps/BHS2018.h
        bdps/BKPX2015.h
        bdps/BSX2009.h
        bdps/Degree3.h
        bdps/KPT2017.h
        bdps/KPX2010.h
        bdps/KX2012.h
        bdps/LW2004.h
        bdps/types.h

        BoundedDegreePlaneSpanners.h

        delaunay/DelaunayL2.h
        delaunay/DelaunayLinf.h
        delaunay/DelaunayTD.h

        constants.h
        types.h

        main.cpp Utilities.h)

target_link_libraries(libspanner gmp)
target_link_libraries(libspanner mpfr)
target_link_libraries(libspanner boost_system)
target_link_libraries(libspanner boost_thread)
#target_link_libraries(libspanner stdc++fs)