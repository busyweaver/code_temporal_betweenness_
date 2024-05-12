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
  std::string cost = "";
  std::string type = "";
  std::string percentTime;
  bool runNonStrict = false;
  bool runBoost = false;
  bool runBFS = false;
  bool runNoLoop = false;
  bool runBellman = false;
  bool edgesDirected = false;
  bool originalNodeIds = false;
  bool readFromFile = false;
  bool runGeneral = true;
  bool runTest = false;
  bool write; // write the results into files
  std::string filename;
  std::string epsilon;
  std::string optimal_cost;
};

// write execution times into a file
void writeTime(const akt::Graph& g, BenchmarkSettings& bs, std::map<std::string,double> ti, std::string s)
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
  file.open (path+"/info"+"_"+s+".txt");
  file <<  "number nodes " << g.N() << "\n";
  file <<  "number events " << g.T() << "\n";
  if (bs.edgesDirected)
    file <<  "number edges " << g.M() << "\n";
  else
    file <<  "number edges " << g.M()/2 << "\n";
  file <<  "number edges_static " << g.MS() << "\n";
  for(auto &e : ti)
      file << e.first << " " << e.second << "\n";
  file.close();
}

// write node ids into a file to recover the nodes to which betweenness values belong
void writeNodeIds(std::vector<std::string> ids, BenchmarkSettings& bs)
{
  std::string s = "nodesIds";
  std::string path;
  if(bs.edgesDirected == false)
    path = bs.filename+"_undirected_exp";
  else
    path = bs.filename+"_directed_exp";
  if(bs.runBoost)
    path = path+"_boost";
  const char* str = path.c_str();

  mkdir(str,0777);

  std::ofstream nodeIds;
  nodeIds.open (path+"/"+ s +".txt");
  for(int i = 0;i< ids.size(); i++)
    nodeIds << i << " "<< ids[i] << "\n";
  nodeIds.close();

}

// write real time stamps to a file in order to recover time betweenness values
void writeTimestamps(std::vector<int> ti, BenchmarkSettings& bs)
{
  std::string s = "timeStamps";
  std::string path;
  if(bs.edgesDirected == false)
    path = bs.filename+"_undirected_exp";
  else
    path = bs.filename+"_directed_exp";
  if(bs.runBoost)
    path = path+"_boost";
  const char* str = path.c_str();

  mkdir(str,0777);

  std::ofstream time;
  time.open (path+"/"+ s +".txt");
  for(int i = 0;i< ti.size(); i++)
    time << i << " "<< ti[i] << "\n";
  time.close();

}

// write static betweenness centrality into a file
void writeStaticBet(std::vector<double> p, BenchmarkSettings& bs)
{
  std::string s = "staticBet";
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

  std::ofstream staticBet;
  staticBet.open (path+"/"+ s +".txt");
  for(auto &e : p)
    staticBet << e << "\n";
  staticBet.close();

}

// write betweenness values into a file
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

// check the values of B(v) with the algorithm of Buss et al. on shortest and shortest foremost paths
bool check_all(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &y)
{
  bool res = true;
  for(int i = 0; i < x.size() ; i++)
    {
      for(int j = 0; j < x[i].size() ; j++)
        { 
      if( abs(x[i][j] - y[i][j]) > 0.01)
        res  = false;
        }
    }

  if(!res)
    {
      for(int i = 0; i < x.size() ; i ++)
        {
          for(int j = 0; j < x[i].size() ; j++)
            { 
              std::cout<< "i,j " << i << " " << j << " "  << x[i][j] << " " << y[i][j] << "\n";
            }
        }
    }
  return res;
}

bool check(std::vector<double> &x, std::vector<double> &y)
{
  bool res = true;
  for(int i = 0; i < x.size() ; i ++)
    if( abs(x[i] - y[i]) > 0.01)
      res  = false;
  if(!res)
    {
      for(int i = 0; i < x.size() ; i ++)
        {
          std::cout << x[i] << " " << y[i] << "\n";
        }
    }
  return res;
}

// runs benchmark on a specific temporal graph with shortest paths variants using a temporal BFS
BenchmarkResults runBenchmarksShort(const akt::Graph& g, BenchmarkSettings& bs)
{

  BenchmarkResults res;
  std::map<std::string,double> ti;
  std::vector<std::pair<std::string,std::string>> cost_type;
  if(bs.cost == "")
    cost_type = {{"shortestrestless","passive"},{"shortestrestless","active"}, {"shortest","passive"},{"shortest","active"}, {"shortestforemost","passive"}};
  else
    {
      cost_type = {{bs.cost,bs.type}};
    }


  if (bs.numberNodes.size() == 0)
    bs.numberNodes = "-1";
  if (bs.percentTime.size() == 0)
    bs.percentTime = "-1";
    for (auto &st: cost_type){
            std::vector<std::string> strict;
            if(bs.runStrict == true)
                strict.push_back("strict");
            else
              strict.push_back("non-strict");

            for (auto &stri: strict){
                bool str_bool;
                if (stri == "strict")
                  str_bool = true;
                else
                  str_bool = false;
                //std::clog << "Starting  " << st.first << " " << st.second << std::endl;
                std::cout << "results saved into files ?" << bs.write << " *0 means no, 1 yes* \n";
                std::cout << stri + "_"+st.first+"_"+st.second << "\n";
                auto start = std::chrono::high_resolution_clock::now();
                std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>> , std::vector<double> > y;
                if(bs.runBoost)
                  y = optimalBetweennessBoost(g, str_bool, st.first, st.second, stoi(bs.numberNodes), stoi(bs.percentTime));
                else
                  y = optimalBetweenness(g, str_bool, st.first, "le", st.second, stoi(bs.numberNodes), bs.runBellman, bs.runBFS, bs.runNoLoop);
                auto end = std::chrono::high_resolution_clock::now();
                std::pair< std::vector<std::vector<double>>, std::vector<std::vector<double>> > x;
                std::vector<std::vector<double>> bet;
                std::vector<std::vector<double>> bet_ex;
                std::vector<double> bet_sum;
                if(bs.write)
                  {
                    std::tie (bet, bet_ex, bet_sum) = y;
                    x.first = bet;
                    x.second = bet_ex;
                    writeToFile(g, bs, stri+"_"+st.first+"_"+st.second+"_"+bs.numberNodes+"_"+bs.percentTime, x); 
                  }
                std::cout << "end calcul "<< "\n" << std::flush;
                std::chrono::duration<double> time = end - start;
                res.optimalTime[stri+"_"+st.first+"_"+st.second] = time.count();
                std::ofstream file;
                if (bs.runBoost)
                  file.open (bs.filename+"_"+bs.cost+"_"+bs.type+"_bfs");
                else
                  {
                    if(bs.runBellman)
                      file.open (bs.filename+"_"+bs.cost+"_"+bs.type+"_bell");
                    else
                      file.open (bs.filename+"_"+bs.cost+"_"+bs.type+"_dij");
                  }
                file << time.count();
                file.close();
                std::cout << "time elapsed : "  << time.count() << "\n";

                // run test if needed, to automatically check values, only ran if program with option -t
                if(bs.runTest && ((st.first == "shortest" && st.second== "passive")  || (st.first == "shortestforemost" && st.second== "passive")) ){
                  std::cout << "START TEST \n";


                  y = optimalBetweenness(g, str_bool, st.first, "le", st.second, stoi(bs.numberNodes), false, false, bs.runNoLoop);
                  std::pair< std::vector<std::vector<double>>, std::vector<std::vector<double>> > x;
                  std::vector<std::vector<double>> bet;
                  std::vector<std::vector<double>> bet2;
                  std::vector<std::vector<double>> bet_ex;
                  std::vector<double> bet_sum;
                  std::tie (bet, bet_ex, bet_sum) = y;

                  auto start = std::chrono::high_resolution_clock::now();
                  auto p = shortestBetweenness(g, str_bool);
                  auto end = std::chrono::high_resolution_clock::now();
                  std::chrono::duration<double> time = end - start;
                  res.optimalTime["buss_"+st.first+"_"+st.second] = time.count();
                  std::cout << "time elapsed Test: "  << time.count() << "\n";
                  bool test;
                  if (st.first == "shortest" && st.second== "passive")
                    test = check(p.first, bet_sum);
                  else
                    test = check(p.second, bet_sum);


                  y = optimalBetweenness(g, str_bool, st.first, "le", st.second, stoi(bs.numberNodes), true, false, bs.runNoLoop);
                  std::tie (bet2, bet_ex, bet_sum) = y;

                  test = check_all(bet,bet2);

                  y = optimalBetweenness(g, str_bool, st.first, "le", st.second, stoi(bs.numberNodes), false, true, bs.runNoLoop);
                  std::tie (bet2, bet_ex, bet_sum) = y;

                  test = check_all(bet,bet2);

                  if(!test){
                    std::cout << "probleme test";
                    //exit(-1);
                  }
                  else
                    std::cout << "test OK! \n";
                }
                if(bs.write)
                  {
                    start = std::chrono::high_resolution_clock::now();
                    auto p = shortestBetweennessStatic(g);
                    end = std::chrono::high_resolution_clock::now();
                    time = end - start;
                    res.optimalTime["static_betweenness"] = time.count();
                    writeStaticBet(p, bs);
                    writeTime(g, bs, res.optimalTime, stri+"_"+st.first+"_"+st.second+"_"+bs.numberNodes+"_"+bs.percentTime);
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
    BenchmarkResults br;
    if(bs.write)
      {
        writeNodeIds(ids,bs);
        writeTimestamps(g.events,bs); 
      }
    br = runBenchmarksShort(g, bs);
    br.inputIds = std::move(ids);
    br.n = g.N();

    return br;
}

void outputBenchmarkResults(const BenchmarkSettings& bs, const BenchmarkResults& br)
{
    std::clog << "Time for algorithms (in seconds):\n";
    for(auto &ot: br.optimalTime)
      std::cout << ot.first << " " << ot.second;
    std::cout <<  "\n";
}

int main (int argc, char** argv)
{


  BenchmarkSettings bs;
  po::options_description desc("usage: btwBenchmark [options]\nRuns the different betweenness centrality algorithms on the graph input via stdin (or file if -f used). Available options");
  desc.add_options()
		("help,h", "write help message")
    ("test,t", "run test on passive shortest and shortest foremost")
		("filename,f", po::value<std::string>(&(bs.filename)), "instead of reading the graph from stdin, use the file given in the argument")
    ("epsilon,e", po::value<std::string>(&(bs.epsilon)), "value for event time difference")
    ("numbernodes,a", po::value<std::string>(&(bs.numberNodes)), "number of nodes treated for heuristic")
    ("percenttime,p", po::value<std::string>(&(bs.percentTime)), "percentage of time looking an integer")
    ("optimal,opt", po::value<std::string>(&(bs.optimal_cost)), "choose a cost function for the general model if not specified shortest paths are selected")
    ("graph-directed,d", "interpret the edges in the graph as directed edges")
    ("write,w", "writes the results into files to save values of betweenness")
		("strict,s", "run the strict versions betweenness algorithm, default non-strict (for BFS version) in other cases just implement the function as strict")
    ("noloop,o", "if predecessor never makes loop, activate to speed up computation")
    ("bellman-ford,z", "run bellman-ford variant betweenness algorithm")
    ("boost,b", "run the accelerated version of shortest")
    ("bfs,x", "run the generic BFS")
    ("cost,c", po::value<std::string>(&(bs.cost)),  "Cost function")
    ("type,y",po::value<std::string>(&(bs.type)), "type active/passive")
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
  bs.runBellman = vm.count("bellman-ford") > 0;
  bs.runNoLoop = vm.count("noloop") > 0;
  bs.write = vm.count("write") > 0;
  bs.runBoost = vm.count("boost") > 0;
  bs.runBFS = vm.count("bfs") > 0;
  bs.runTest = vm.count("test") > 0;
  bs.originalNodeIds = vm.count("originalNodeIds") > 0;
    try {
        auto br = readGraphRunBenchmarks(bs);
        outputBenchmarkResults(bs, br);
    } catch (...) {
        return 1;
    }
}
