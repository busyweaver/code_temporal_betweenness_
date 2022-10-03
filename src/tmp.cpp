

#include<networkit/graph/Graph.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <set>
#include <networkit/distance/BFS.hpp>

int main(int argc, char* argv[]) {

  // create a graph Graph(count n = 0, bool weighted = false, bool directed = false, bool edgesIndexed = false)
  NetworKit::Graph G = NetworKit::Graph(0, false, true, false);
  NetworKit::node no = G.addNode();
  NetworKit::node no2 = G.addNode();
  NetworKit::node no3 = G.addNode();
  //networkit::Node no2 = G.addNode();
  G.addEdge(no,no2);
  G.addEdge(no2,no);
  G.addEdge(no,no3);
  G.forNodes([&](NetworKit::node u) { printf("node %ld\n",u); });
  NetworKit::StronglyConnectedComponents scc(G);
  scc.run();
  const auto partition = scc.getPartition();
  const auto par = scc.getPartition().getSubsets();
  std::set<std::set<NetworKit::index>>::iterator itr;
  const auto sizes = scc.getPartition().subsetSizes();
  NetworKit::Graph GG = NetworKit::Graph(sizes.size(), false, true, false);
  // Displaying set elements
  for (itr = par.begin(); itr != par.end(); itr++)
    {
          std::set<NetworKit::index>::iterator itr2;
          itr2 = (*itr).begin();
          printf("from %ld\n", *itr2);
          NetworKit::BFS bfs(G, (*itr2), false, true);
          bfs.run();
          std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
          //          std::vector<NetworKit::node> vi;
          for(NetworKit::node i : stack)
            {
              if( partition.subsetOf(i) !=  partition.subsetOf(*itr2))
                {
                  GG.addEdge(partition.subsetOf(*itr2),partition.subsetOf(i));
                }
            }
          //std::cout << *itr2 << " ";
    }
  printf("Edges of GG\n");
  GG.forEdges([&](NetworKit::node u, NetworKit::node v) { printf("edge %ld,%ld\n",u, v); });
  printf("sub of %ld\n", partition.subsetOf(0));
  //std::set<std::set<index>> 

  const auto part_vec = scc.getPartition().getVector(); 
  const auto ids = scc.getPartition().getSubsetIds();
  for (auto const &i: sizes) {
    std::cout << i << std::endl;
  }
  printf("part_vector");
  for (auto const &i: part_vec) {
    std::cout << i << std::endl;
  }
  printf("ids");
  for (auto const &i: ids) {
    std::cout << i << std::endl;
  }


}
