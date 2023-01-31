#pragma once
#include "paths.h"
#include <unordered_set>
#include<iostream>
struct VertexAppearance
{
  // Vertex id
  int v;
  // Timestamp
  int time;
  bool operator<(const VertexAppearance & other) const
  {
    std::pair<int,int> p2(other.v,other.time);
    std::pair<int,int> p (v,time);

    return  p < p2;
  }
};

// Stores data about a vertex appearance
bool operator==(const VertexAppearance& lhs, const VertexAppearance& rhs);
// Simple hashing function for akt::VertexAppearance (so that we can use it with std::unordered_set)
namespace std
{
  template<> struct hash<VertexAppearance>
  {
    std::size_t operator()(const VertexAppearance& va) const noexcept
    {
      auto h1 = std::hash<int>{}(va.v);
      auto h2 = std::hash<int>{}(va.time);
      return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
    }
  };
}

// Adapted from https://en.cppreference.com/w/cpp/utility/hash
  struct OptimalBetweennessData
  {
  OptimalBetweennessData(int n, int T)
    : deltasvvt{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
      deltadot{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
      betweenness{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
      betweenness_exact{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
      sigma{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0)) },
      sigmadot{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0)) },
      totalSigmaT{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
      pre{ std::vector<std::vector<std::unordered_set<VertexAppearance>>>(n, std::vector<std::unordered_set<VertexAppearance>>(T)) },
       cur_best{ std::vector<std::vector<double>>(n, std::vector<double>(T, std::numeric_limits<double>::infinity())) },
      opt_walk{ std::vector<std::vector<Path*>>(n, std::vector<Path*>(T, (Path*) nullptr) ) },
      optimalNode{ std::vector<double>(n, std::numeric_limits<double>::infinity()) },
      totalSigma{ std::vector<double>(n, 0) },
      totalBetweenness{ std::vector<double>(n, 0.0) }
  { }
  OptimalBetweennessData(const akt::Graph& g)
    : OptimalBetweennessData(g.N(), g.T())
  { }

  // The \delta_{s\cdot} array on each iteration
  std::vector<std::vector<double>> deltadot;
  std::vector<std::vector<double>> betweenness;
    std::vector<std::vector<double>> betweenness_exact;
  std::vector<std::vector<double>> deltasvvt;
  // Numbers of shortest paths from a source to each vertex appearance
  std::vector<std::vector<double>> sigma;
  std::vector<std::vector<double>> sigmadot;
  std::vector<std::vector<int>> totalSigmaT;
  // Sets of predecessors for each vertex appearance
  std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre;
  // Distances to each vertex appearance
  std::vector<std::vector<double>> cur_best;
  std::vector<std::vector<Path*>> opt_walk;
  // Distances to each vertex (*not* vertex __appearance__)
  std::vector<double> optimalNode;
  // Number of shortests paths from a source to each vertex (*not* vertex __appearance__)
  std::vector<double> totalSigma;
    std::vector<double> totalBetweenness;
};
