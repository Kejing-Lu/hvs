cmake_minimum_required (VERSION 2.6)
project (hnsw_lib)

find_package(OpenCV REQUIRED)

include_directories("${PROJECT_BINARY_DIR}" ${OpenCV_INCLUDE_DIRS})



set(SOURCE_EXE main.cpp)            

set(SOURCE_LIB sift_1b.cpp)

add_library(sift_test STATIC ${SOURCE_LIB})


add_executable(main ${SOURCE_EXE})
SET( CMAKE_CXX_FLAGS  "-Ofast -lrt -DNDEBUG -std=c++11 -I /usr/local/boost_1_75_0 -DHAVE_CXX0X -openmp -march=native -fpic -w -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=0" )
target_link_libraries(main sift_test ${OpenCV_LIBS})  
