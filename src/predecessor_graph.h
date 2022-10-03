#pragma once

#include<networkit/graph/Graph.hpp>
#include "algorithms.h"

class Predecessor {
 public:
  NetworKit::Graph& g;
  std::unordered_set<int> sources;
  std::unordered_set<int> sinks;
  std::map<int, std::map<int, std::vector<int> > > ordered_neighb;
  Predecessor () {}


};


Predecessor* PredecessorGraph(const akt::Graph& g, std::vector<std::map<int, std::unordered_set<akt::VertexAppearance>>> pre, int node);
std::pair<std::unordered_set<int>, std::unordered_set<int>> RemoveInfiniteFromPredecessor(int s, Predecessor& G, OptimalBetweennessData& sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), std::string walk_type, std::vector<int> &events, std::vector<int>& events_rev);
unordered_set<akt::VertexAppearance> Infinite_closure(Predecessor& G, std::vector<int> &events, std::map<int, int> &events_rev, OptimalBetweennessData &sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), double n, std::unordered_set<int> &node_inf);
void VolumePathAt(Predecessor& G, int s, set<NetworKit::index> sinks, OptimalBetweennessData &sbd);
void OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& s, double n, double (*cost)(Path, int, double))
void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd);
void PredecessorGraphToOrdered(Predecessor& G, int  T, std::map<int, int> ev_rev);
