#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <tuple>

#include <boost/program_options.hpp>
#include "algorithms.h"

namespace po = boost::program_options;

struct BenchmarkResults 
{
    std::vector<double> shortest, foremost, strictShortest, strictForemost, prefix;
  std::vector<double> optimal;
    int n;
    std::vector<std::string> inputIds;
    double nonStrictTime = -1.0;
    double strictTime = -1.0;
    double prefixTime = -1.0;
  double optimalTime = -1.0;
};

struct BenchmarkSettings
{
    bool runStrict = true;
    bool runNonStrict = true;
    bool runPrefix = true;
    bool edgesDirected = false;
    bool originalNodeIds = false;
    bool readFromFile = false;
  bool runGeneral = true;
    std::string filename;
};

BenchmarkResults runBenchmarks(const akt::Graph& g, const BenchmarkSettings& bs)
{
    BenchmarkResults res;

    if (bs.runNonStrict) {
        std::clog << "Starting non-strict shortest / shortest foremost." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::tie(res.shortest, res.foremost) = shortestBetweenness(g, false);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        res.nonStrictTime = time.count();
    }
    if (bs.runStrict) {
        std::clog << "Starting strict shortest / shortest foremost." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::tie(res.strictShortest, res.strictForemost) = shortestBetweenness(g, true);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        res.strictTime = time.count();
    }
    if (bs.runPrefix) {
        std::clog << "Starting prefix foremost." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        res.prefix = prefixForemostBetweenness(g);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        res.prefixTime = time.count();
    }

    if (bs.runGeneral) {
      std::clog << "Starting general." << std::endl;
      auto start = std::chrono::high_resolution_clock::now();
      res.optimal = optimalBetweenness(g, false, "shortestfastest", "le", "passive");
      printf("finiiiiiiiiiii\n");
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time = end - start;
      res.optimalTime = time.count();
    }

    return res;
}

auto readGraphFromFile(const BenchmarkSettings& bs)
{
    std::ifstream ifs{ bs.filename };
    if (!ifs) {
        std::cout << "Error trying to open the file \"" << bs.filename << "\"\n";
        throw 1;
    }
    return akt::readReduceGraph(ifs, bs.edgesDirected);
}

BenchmarkResults readGraphRunBenchmarks(const BenchmarkSettings& bs)
{
    // Not the perfect solution from the perspective of error handling (opening the file and not checking for success), but simple to write cleanly)
    auto [g, ids] = bs.readFromFile ? readGraphFromFile(bs) : akt::readReduceGraph(std::cin, bs.edgesDirected);
    std::clog << "Graph read: " << g.N() << " nodes, ";
    if (bs.edgesDirected)
        std::clog << g.M() << " (unique) directed edges, ";
    else
        std::clog << g.M() / 2 << " (unique) edges, ";
    std::clog << g.T() << " (non-empty) timesteps\n";
    auto br = runBenchmarks(g, bs);
    br.inputIds = std::move(ids);
    br.n = g.N();

    return br;
}

void outputBenchmarkResults(const BenchmarkSettings& bs, const BenchmarkResults& br)
{
    std::clog << "Time for algorithms (in seconds):\n";
    std::clog << "Non-strict, strict, prefix foremost, optimal general time\n";
    std::clog << br.nonStrictTime << ", " << br.strictTime << ", " << br.prefixTime << br.optimalTime<< '\n';
    std::cout << br.nonStrictTime << ", " << br.strictTime << ", " << br.prefixTime <<" " <<'\n';
    std::clog << "Computed betweenness measures:\n";
    std::clog << "Node, non-strict shortest, non-strict shortest foremost, strict shortest, strict shortest foremost, optimal general, prefix betweenness\n";
    for (int i = 0; i < br.n; ++i) {
        if (bs.originalNodeIds)
            std::cout << br.inputIds[i];
        else
            std::cout << i;
        std::cout << ", ";
        if (bs.runNonStrict)
            std::cout << br.shortest[i] << ", " << br.foremost[i] << ", ";
        else
            std::cout << "-1, -1, ";
        if (bs.runStrict)
            std::cout << br.strictShortest[i] << ", " << br.strictForemost[i] << ", ";
        else
            std::cout << "-1, -1, ";
        if (bs.runGeneral)
          std::cout << br.optimal[i] << ", ";
        else
          std::cout << "-1\n";
        if (bs.runPrefix)
            std::cout << br.prefix[i] << '\n';
        else
            std::cout << "-1\n";
    }
}

int main (int argc, char** argv)
{
    BenchmarkSettings bs;
    po::options_description desc("usage: btwBenchmark [options]\nRuns the different betweenness centrality algorithms on the graph input via stdin (or file if -f used). Available options");
    desc.add_options()
		("help,h", "write help message")
		("filename,f", po::value<std::string>(&(bs.filename)), "instead of reading the graph from stdin, use the file given in the argument")
        ("graph-directed,d", "interpret the edges in the graph as directed edges")
		("no-strict,s", "don't run the strict shortest (foremost) betweenness algorithm")
		("no-non-strict,n", "don't run the non-strict shortest (foremost) betweenness algorithm")
		("no-prefix,p", "don't run the the strict prefix-foremost betweenness algorithm")
		("original-node-ids,u", "in the output, use original node ids instead of the arbitrary integer ids assigned by the program");
        
    po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << '\n';
        return 0;
    }
    bs.readFromFile = vm.count("filename") > 0;
    bs.edgesDirected = vm.count("graph-directed") > 0;
    bs.runStrict = vm.count("no-strict") <= 0;
    bs.runNonStrict = vm.count("no-non-strict") <= 0;
    bs.runPrefix = vm.count("no-prefix") <= 0;
    bs.originalNodeIds = vm.count("originalNodeIds") > 0;
    
    try {
        auto br = readGraphRunBenchmarks(bs);
        outputBenchmarkResults(bs, br);
    } catch (...) {
        return 1;
    }
}
