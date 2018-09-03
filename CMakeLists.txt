cmake_minimum_required( VERSION 2.8 )
project ( GPS )

set( CMAKE_BUILD_TYPE "Debug" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )


list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

############### dependencies ######################
#boost
find_package(Boost REQUIRED COMPONENTS system)
include_directories( ${Boost_INCLUDE_DIR} )
#threads
find_package(Threads)

# Eigen
include_directories( "/usr/local/include/eigen3" )
# OpenCV
#find_package( OpenCV 3.4 REQUIRED )
#include_directories( ${OpenCV_INCLUDE_DIRS} )

#set( THIRD_PARTY_LIBS ${OpenCV_LIBS} )
############### dependencies ######################

set(SOURCE_FILES
        SVs.cpp
        NtripRTK.cpp
        PosSolver.cpp
        SerialData.cpp
        GNSS.cpp
        GnssTime.cpp GnssTime.h EphemSp3.cpp EphemSp3.h CommonInclude.h CommonInclude.cpp)
add_executable(GPS  main.cpp ${SOURCE_FILES} CommonInclude.h)
add_executable(testrtk testrtk.cpp ${SOURCE_FILES} )
#target_link_libraries( SYS_GPS ${OpenCV_LIBS} )
target_link_libraries( GPS  ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#target_link_libraries(SYS_GPS )

target_link_libraries(testrtk ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
