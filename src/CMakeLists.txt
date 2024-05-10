cmake_minimum_required (VERSION 2.8)
project (temporal-betweenness)
include (cmake/libs.cmake)

set(CMAKE_CXX_FLAGS " -std=c++17 -pedantic -O3 ${CMAKE_CXX_FLAGS}")

#find_package(Boost 1.70.0)
#if (Boost_FOUND)
#	include_directories(${Boost_INCLUDE_DIR}/boost)
#else ()
install_boost(gcc program_options)
#endif ()

add_executable(btwBenchmark
	btwBenchmark.cpp
	algorithms2.cpp
	graph.cpp
  paths.cpp
  predecessor_graph.cpp
  vertex_appearance.cpp
)
# if networkit is installed on globaly put this
target_link_libraries(btwBenchmark boost_program_options -lnetworkit -g)
#target_link_libraries(btwBenchmark boost_program_options -L/home/naima/include -lnetworkit -g)

