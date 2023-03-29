#pragma once
#include "paths.h"
#include "predecessor_graph.h"
#include <unordered_set>

namespace akt {
	// Returns the shortest betweenness and shortest foremost betweenness (with strictness depending on the parameter)
	std::pair<std::vector<double>, std::vector<double>> shortestBetweenness(const akt::Graph& g, bool strict);
  std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>> , std::vector<double> > optimalBetweenness(const Graph& g, bool strict, std::string cost, std::string cmp, std::string walk_type, int numberNodes);
  std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>> , std::vector<double> > optimalBetweennessBoost(const Graph& g, bool strict, std::string rest, std::string walk_type, int numberNodes);

	// Returns the strict prefix-foremost betwenness
	std::vector<double>	prefixForemostBetweenness(const akt::Graph& g);

	// Returns the shortest betweenness for a static graph (i.e. all edges have timestamp 0) using Brandes' algorithm
	std::vector<double> shortestBetweennessStatic(const akt::Graph& g);


} // end namespace akt

  // struct VertexAppearance
  // {
  //   // Vertex id
  //   int v;
  //   // Timestamp
  //   int time;
  //   bool operator<(const VertexAppearance & other) const
  //   {
  //     std::pair<int,int> p2(other.v,other.time);
  //     std::pair<int,int> p (v,time);

  //     return  p < p2;
  //   }
  // };
  // // Stores data about a vertex appearance

  // bool operator==(const VertexAppearance& lhs, const VertexAppearance& rhs)
  // {
  //   return (lhs.v == rhs.v) && (lhs.time == rhs.time);
  // }
  // // Simple hashing function for akt::VertexAppearance (so that we can use it with std::unordered_set)
  // namespace std
  // {
  //   template<> struct hash<VertexAppearance>
  //   {
  //     std::size_t operator()(const VertexAppearance& va) const noexcept
  //     {
  //       auto h1 = std::hash<int>{}(va.v);
  //       auto h2 = std::hash<int>{}(va.time);
  //       return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
  //     }
  //   };
  // }
  // // Adapted from https://en.cppreference.com/w/cpp/utility/hash

//   struct OptimalBetweennessData
//   {
//   OptimalBetweennessData(int n, int T)
//     : deltasvvt{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
//       deltadot{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
//       betweenness{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
//       sigma{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
//       sigmadot{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
//       totalSigmaT{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
//       pre{ std::vector<std::vector<std::unordered_set<VertexAppearance>>>(n, std::vector<std::unordered_set<VertexAppearance>>(T)) },
//       cur_best{ std::vector<std::vector<double>>(n, std::vector<double>(T, std::numeric_limits<double>::infinity())) },
//       opt_walk{ std::vector<std::vector<Path>>(n, std::vector<Path>(T, Path(NULL,NULL))) },
//       optimalNode{ std::vector<double>(n, std::numeric_limits<double>::infinity()) },
//       totalSigma{ std::vector<int>(n, 0) },
//       stack{ std::vector<VertexAppearance>() }
//   { }
//   OptimalBetweennessData(const akt::Graph& g)
//     : OptimalBetweennessData(g.N(), g.T())
//   { }

//   // The \delta_{s\cdot} array on each iteration
//   std::vector<std::vector<double>> deltadot;
//   std::vector<std::vector<double>> betweenness; 
//   std::vector<std::vector<double>> deltasvvt;
//   // Numbers of shortest paths from a source to each vertex appearance
//   std::vector<std::vector<int>> sigma;
//   std::vector<std::vector<int>> sigmadot;
//   std::vector<std::vector<int>> totalSigmaT;
//   // Sets of predecessors for each vertex appearance
//   std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre;
//   // Distances to each vertex appearance
//   std::vector<std::vector<double>> cur_best;
//   std::vector<std::vector<Path>> opt_walk;
//   // Distances to each vertex (*not* vertex __appearance__)
//   std::vector<double> optimalNode;
//   // Number of shortests paths from a source to each vertex (*not* vertex __appearance__)
//   std::vector<int> totalSigma;
//   // Stack of vertex apperances in order of their discovery with the bfs
//   std::vector<VertexAppearance> stack;
// };
