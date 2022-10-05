cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
include(ExternalProject)

macro(install_boost boost_toolset boost_libs)
ExternalProject_Add(boost
	URL https://boostorg.jfrog.io/artifactory/main/release/1.70.0/source/boost_1_70_0.tar.gz
	URL_HASH SHA256=882b48708d211a5f48e60b0124cf5863c1534cd544ecd0664bb534a4b5d506e9
	CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}  ./bootstrap.sh --with-toolset=${boost_toolset} --prefix=${CMAKE_CURRENT_BINARY_DIR}/boost --with-libraries=${boost_libs}
	BUILD_COMMAND ./b2  variant=release link=static runtime-link=static threading=multi cxxflags=${CMAKE_CXX_FLAGS}
	INSTALL_COMMAND ./b2 variant=release link=static runtime-link=static threading=multi cxxflags=${CMAKE_CXX_FLAGS} install 	
	BUILD_IN_SOURCE 1
)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/boost/include)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/boost/lib)
endmacro(install_boost)





