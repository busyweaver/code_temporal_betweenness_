#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <tuple>
#include <iostream>
//#include <filesystem>

#include <vector>
#include <string>

#include <boost/program_options.hpp>
#include "algorithms.h"
#include <sys/stat.h>


namespace po = boost::program_options;

struct BenchmarkResults 
{
  int n;
  std::map<std::string, std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>>> results_general;
  std::vector<std::string> inputIds;
  std::map<std::string,double> optimalTime;
};

struct BenchmarkSettings
{
  bool runStrict = false;
  bool edgesDirected = false;
  bool originalNodeIds = false;
  bool readFromFile = false;
  bool is_epsilon = false;
  bool runGeneral = true;
  std::string filename;
  std::string epsilon;
  std::string optimal_cost;
};

void writeToFile(const akt::Graph& g, BenchmarkSettings& bs, std::string s, std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> p)
{
  std::string path = bs.filename+"_exp";
  const char* str = path.c_str();

  mkdir(str,0777);

  std::ofstream file_bet;
  std::ofstream file_bet_exact;
  file_bet.open (path+"/bet_"+s+".txt");
  file_bet_exact.open (path+"/bet_exact_"+s+".txt");
  for(int i = 0; i < g.N(); i++)
    {
      for(int j = 0; j < g.T(); j ++)
        {
          file_bet << p.first[i][j] << "\n";
          file_bet_exact << p.second[i][j] << "\n";
        }
    }
  file_bet.close();
  file_bet_exact.close();

}

BenchmarkResults runBenchmarks(const akt::Graph& g, BenchmarkSettings& bs)
{
    BenchmarkResults res;

    std::vector<std::pair<std::string,std::string>> cost_type{{"shortest","passive"}, {"shortest","active"}, {"shortestfastest","passive"} , {"shortestfastest","active"}, {"foremost","passive"} , {"shortestforemost","passive"}};
    //std::vector<std::pair<std::string,std::string>> cost_type{{"shortest","passive"}};
    for (auto &st: cost_type)
      {
            std::vector<std::string> strict;
            strict.push_back("false");
            if(bs.runStrict == true)
                strict.push_back("true");


            for (auto &stri: strict)
              {
                std::clog << "Starting  " << st.first << " " << st.second << std::endl;
                auto start = std::chrono::high_resolution_clock::now();
                if(stri == "false")
                  {
                    std::cout << "non-strict_"+st.first+"_"+st.second;
                    //                   res.results_general["non-strict_"+st.first+"_"+st.second] = optimalBetweenness(g, false, st.first, "le", st.second);
                    writeToFile(g, bs, "non-strict_"+st.first+"_"+st.second, optimalBetweenness(g, false, st.first, "le", st.second));
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["non-strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << " "  << time.count() << "\n";
                  }
                else
                  {
                    std::cout << "strict_"+st.first+"_"+st.second;
                    //res.results_general["strict_"+st.first+"_"+st.second] = optimalBetweenness(g, true, st.first, "le", st.second);
                    writeToFile(g, bs, "non-strict_"+st.first+"_"+st.second, optimalBetweenness(g, true, st.first, "le", st.second));
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << " "  << time.count() << "\n";
                  }

              }
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
    std::string epsilon;
    if (bs.epsilon.size() == 0)
      epsilon = "0";
    else
      epsilon = bs.epsilon;
    std::cout << epsilon << "<- \n";
    return akt::readReduceGraph(ifs, bs.edgesDirected, stod(epsilon));
    //return akt::readReduceGraph(ifs, bs.edgesDirected);
}

BenchmarkResults readGraphRunBenchmarks(BenchmarkSettings& bs)
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
    for(auto &ot: br.optimalTime)
      std::cout << ot.first << " " << ot.second;
}

int main (int argc, char** argv)
{
  BenchmarkSettings bs;
  po::options_description desc("usage: btwBenchmark [options]\nRuns the different betweenness centrality algorithms on the graph input via stdin (or file if -f used). Available options");
  desc.add_options()
		("help,h", "write help message")
		("filename,f", po::value<std::string>(&(bs.filename)), "instead of reading the graph from stdin, use the file given in the argument")
    ("epsilon,e", po::value<std::string>(&(bs.epsilon)), "value for event time difference")
    ("optimal,opt", po::value<std::string>(&(bs.optimal_cost)), "choose a cost function for the general model if not specified shortest paths are selected")
    ("graph-directed,d", "interpret the edges in the graph as directed edges")
		("strict,s", "run the strict versions betweenness algorithm")
		("no-non-strict,n", "don't run the non-strict versions betweenness algorithm");
  po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 0;
  }
  bs.readFromFile = vm.count("filename") > 0;
  bs.edgesDirected = vm.count("graph-directed") > 0;
  bs.is_epsilon = vm.count("epsilon") > 0;
  bs.runStrict = vm.count("strict") > 0;
  bs.originalNodeIds = vm.count("originalNodeIds") > 0;
    try {
        auto br = readGraphRunBenchmarks(bs);
        outputBenchmarkResults(bs, br);
    } catch (...) {
        return 1;
    }
}
