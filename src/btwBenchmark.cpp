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
  std::string numberNodes;
  bool runNonStrict = false;
  bool runBoost = false;
  bool edgesDirected = false;
  bool originalNodeIds = false;
  bool readFromFile = false;
  bool runGeneral = true;
  std::string filename;
  std::string epsilon;
  std::string optimal_cost;
};

void writeTime(const akt::Graph& g, BenchmarkSettings& bs, std::map<std::string,double> ti)
{
  std::cout << "write time start "<< "\n" << std::flush;
  std::string path;
  if(bs.edgesDirected == false)
    path = bs.filename+"_undirected_exp";
  else
    path = bs.filename+"_directed_exp";
  if(bs.runBoost)
    path = path+"_boost";
  const char* str = path.c_str();
  mkdir(str,0777);
  std::ofstream file;
  file.open (path+"/info"+"_"+bs.numberNodes+".txt");
  file <<  "number nodes " << g.N() << "\n";
  file <<  "number events " << g.T() << "\n";
  if (bs.edgesDirected)
    file <<  "number edges " << g.M() << "\n";
  else
    file <<  "number edges " << g.M()/2 << "\n";
  for(auto &e : ti)
      file << e.first << " " << e.second << "\n";
  file.close();
}

void writeToFile(const akt::Graph& g, BenchmarkSettings& bs, std::string s, std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> p)
{
  std::cout << "write start "<< "\n" << std::flush;
  std::string path;
  if(bs.edgesDirected == false)
    path = bs.filename+"_undirected_exp";
  else
    path = bs.filename+"_directed_exp";
  if(bs.runBoost)
    path = path+"_boost";
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

BenchmarkResults runBenchmarksShort(const akt::Graph& g, BenchmarkSettings& bs)
{

  BenchmarkResults res;
  std::map<std::string,double> ti;
  std::vector<std::pair<std::string,std::string>> cost_type{{"shortest","passive"}, {"shortest","active"}};
  if (bs.numberNodes.size() == 0)
    bs.numberNodes = "-1";
    for (auto &st: cost_type)
      {
            std::vector<std::string> strict;
            if(bs.runNonStrict == true)
              strict.push_back("false");
            if(bs.runStrict == true)
                strict.push_back("true");


            for (auto &stri: strict)
              {
                //std::clog << "Starting  " << st.first << " " << st.second << std::endl;

                if(stri == "false")
                  {
                    std::cout << "non-strict_"+st.first+"_"+st.second << "\n";
                    auto start = std::chrono::high_resolution_clock::now();
                    auto x = optimalBetweennessBoost(g, false, st.first, st.second, stoi(bs.numberNodes));
                    auto end = std::chrono::high_resolution_clock::now();
                    writeToFile(g, bs, "non-strict_"+st.first+"_"+st.second+"_"+bs.numberNodes, x);
                    std::cout << "end calcul "<< "\n" << std::flush;
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["non-strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << "time elapsed : "  << time.count() << "\n";
                  }
                else
                  {
                    std::cout << "strict_"+st.first+"_"+st.second << "\n";
                    auto start = std::chrono::high_resolution_clock::now();
                    auto x = optimalBetweennessBoost(g, true, st.first, st.second, stoi(bs.numberNodes));
                    auto end = std::chrono::high_resolution_clock::now();
                    writeToFile(g, bs, "strict_"+st.first+"_"+st.second+"_"+bs.numberNodes, x);
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << " "  << time.count() << "\n";
                  }

              }

      }
    writeTime(g, bs, res.optimalTime);
    //if needed uncomment to check
    // std::cout << "verification with shortest passive and  shortest passive"<< "\n" << std::flush;
    // auto start = std::chrono::high_resolution_clock::now();
    // auto p = shortestBetweenness(g, false);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time = end - start;
    // std::cout << "shortest and shortest foremost both passive time  "  << time.count() << "\n";
    // for(auto &elem : p.first)
    //   std::cout << "shortest "<< elem << " "<<"\n";

    return res;
}

BenchmarkResults runBenchmarks(const akt::Graph& g, BenchmarkSettings& bs)
{

  BenchmarkResults res;
  std::map<std::string,double> ti;
  std::vector<std::pair<std::string,std::string>> cost_type{{"shortest","passive"}, {"shortest","active"}, {"shortestfastest","passive"} , {"shortestfastest","active"}, {"foremost","passive"} , {"shortestforemost","passive"}};
  if (bs.numberNodes.size() == 0)
    bs.numberNodes = "-1";
  //  std::vector<std::pair<std::string,std::string>> cost_type{{"shortest","active"}};
    for (auto &st: cost_type)
      {
            std::vector<std::string> strict;
            if(bs.runNonStrict == true)
              strict.push_back("false");
            if(bs.runStrict == true)
                strict.push_back("true");


            for (auto &stri: strict)
              {
                //std::clog << "Starting  " << st.first << " " << st.second << std::endl;

                if(stri == "false")
                  {
                    std::cout << "non-strict_"+st.first+"_"+st.second << "\n";
                    auto start = std::chrono::high_resolution_clock::now();
                    auto x = optimalBetweenness(g, false, st.first, "le", st.second, stoi(bs.numberNodes));
                    auto end = std::chrono::high_resolution_clock::now();
                    writeToFile(g, bs, "non-strict_"+st.first+"_"+st.second+"_"+bs.numberNodes, x);
                    std::cout << "end calcul "<< "\n" << std::flush;
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["non-strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << "time elapsed : "  << time.count() << "\n";
                  }
                else
                  {
                    std::cout << "strict_"+st.first+"_"+st.second << "\n";
                    auto start = std::chrono::high_resolution_clock::now();
                    auto x = optimalBetweenness(g, true, st.first, "le", st.second, stoi(bs.numberNodes));
                    auto end = std::chrono::high_resolution_clock::now();
                    writeToFile(g, bs, "strict_"+st.first+"_"+st.second+"_"+bs.numberNodes, x);
                    std::chrono::duration<double> time = end - start;
                    res.optimalTime["strict_"+st.first+"_"+st.second] = time.count();
                    std::cout << " "  << time.count() << "\n";
                  }

              }

      }
    writeTime(g, bs, res.optimalTime);
    //if needed uncomment to check
    // std::cout << "verification with shortest passive and  shortest foremost passive"<< "\n" << std::flush;
    // auto start = std::chrono::high_resolution_clock::now();
    // auto p = shortestBetweenness(g, false);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time = end - start;
    // std::cout << "shortest and shortest foremost both passive time  "  << time.count() << "\n";
    // for(auto &elem : p.first)
    //   std::cout << "shortest "<< elem << " "<<"\n";
    // for(auto &elem : p.second)
    //   std::cout << "foremost "<< elem << " "<<"\n";

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
    BenchmarkResults br;
    if(bs.runBoost)
      br = runBenchmarksShort(g, bs);
    else
      br = runBenchmarks(g, bs);
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
    ("numbernodes,a", po::value<std::string>(&(bs.numberNodes)), "number of nodes treated for heuristic")
    ("optimal,opt", po::value<std::string>(&(bs.optimal_cost)), "choose a cost function for the general model if not specified shortest paths are selected")
    ("graph-directed,d", "interpret the edges in the graph as directed edges")
		("strict,s", "run the strict versions betweenness algorithm")
    ("boost,b", "run the accelerated version of shortest")
    //    ("all,a", "run the active version on all temporal nodes, if not set run it only on vertex appearances")
		("non-strict,n", "run the non-strict versions betweenness algorithm");
  po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 0;
  }
  bs.readFromFile = vm.count("filename") > 0;
  bs.edgesDirected = vm.count("graph-directed") > 0;
  bs.runStrict = vm.count("strict") > 0;
  bs.runNonStrict = vm.count("non-strict") > 0;
  bs.runBoost = vm.count("boost") > 0;
  bs.originalNodeIds = vm.count("originalNodeIds") > 0;
    try {
        auto br = readGraphRunBenchmarks(bs);
        outputBenchmarkResults(bs, br);
    } catch (...) {
        return 1;
    }
}
