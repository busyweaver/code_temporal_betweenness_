#pragma once

#include <networkit/graph/Graph.hpp>
#include "vertex_appearance.h"
//#include "paths.h"
//#include "graph.h"

class Predecessor {
 public:
  NetworKit::Graph  g;
  std::unordered_set<int> sources;
  std::unordered_set<int> sinks;
  std::map<int, std::vector<int>> times_ord;
  std::map<int, std::map<int, std::vector<int> > > ordered_neighb;

  Predecessor ();
  Predecessor(const akt::Graph& g, std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre, int node);
};



std::pair<std::unordered_set<int>, std::unordered_set<int>> RemoveInfiniteFromPredecessor(int s, Predecessor& G, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph &), bool (*cmp)(double, double), std::string walk_type, const akt::Graph & g);


void VolumePathAt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph &g);

void OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& g, double (*cost)(Path*, int, const akt::Graph&), std::unordered_set<int> node_inf);

void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph&);
void PredecessorGraphToOrdered(Predecessor& G, int  T, int n);


void GeneralContribution(const akt::Graph& g, Predecessor& G, int s, OptimalBetweennessData& sbd , std::map<int, int> &preced, std::string walk_type);

void IntermediaryNodes(int vt, int vtp, std::map<int,int> before,OptimalBetweennessData& sbd, const akt::Graph& g, double s, int pred_time);

std::map<int,int> BeforeNodes(Predecessor& G, const akt::Graph& g);

void CompleteDelta(Predecessor& G, OptimalBetweennessData &sbd, const akt::Graph& g);
