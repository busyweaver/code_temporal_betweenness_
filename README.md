# Temporal Betweenness Algorithms Implementation

## Boost download
Check that cmake version is >= 3.5
replace in src/cmake/libs.cmake the link https://dl.bintray.com/boostorg/release/1.70.0/source/boost_1_70_0.tar.gz link with https://boostorg.jfrog.io/artifactory/main/release/1.70.0/source/boost_1_70_0.tar.gz

## Networkit
### If installed globally
Need to add the path to libnetworkit.so to LD_LIBRARY_PATH
### If installed Locally
Add networkit library to C_INCLUDE_PATH and CPLUS_INCLUDE_PATH I also had to compile (with cmake and then make) ttmath and tlx manually (they are in networkit folder) and add their path C_INCLUDE_PATH and CPLUS_INCLUDE_PATH.
Also need to add the path to libnetworkit.so to LD_LIBRARY_PATH
Finally, in CMakeLists.txt comment the line of options and uncomment the previous one 


## How to build

- `mkdir build.build`
- `cd build.build`
- `cmake ../src`
- `make`

## Structure of the project

The project consists of: 
 - An implementation of a suitable graph representation (`graph.*`)
 - An implementation of temporal betweenness centrality algorithms (`algorithms.*`)
 - Benchmarking program for testing the algorithms on given networks (`btwBenchmark.cpp`)

## How to use

 - Build (see `How to build` above)
 - Use the `./btwBenchmark` executable by either supplying a graph in the proper format (described below) to stdin or in a file by using the -f option
 - Use `./btwBenchmark -h` to see a list of available flags and options



## Graph format for input into the benchmark

Temporal graphs which are read by the benchmark suite need to have the following form:
 - A graph is represented by a sequence of lines, with each corresponding to a (possibly symmetric) temporal edge in the graph
 - Node IDs can be arbitrary, so both `42` as well as `b6f` are valid node IDs. However, note that after reading the graph each node will be assigned an arbitrary integer ID and as such won't be immediately identifiable in the results. To get the results with the original labels use the -u flag
 - Timestamps must be integers, possibly negative, but such that they fit inside 32-bit signed integer
 - Each line of the input must start with the following description of an edge: ID of the origin (tail) node, ID of the destination (head) node, timestamp, all separated by (non-newline) whitespace. However, arbitrary contents may follow and are ignored
 - If the `-d` flag is used, the edges will be considered directed. Otherwise they are assumed to be undirected
 - Self-loops and duplicate edges are allowed, however they will be ignored
Examples of valid graphs:
    0 1 0
    1 2 0
    1 2 1
    2 3 2
Or:
    a 13 -33
    a 13 1
    13 f5 2
    a f5 -1
    zfyk 34$ -0
