//#include "algorithms.h"
#include "fibonacciheap.h"

#include "predecessor_graph.h"
#include <numeric>
#include <queue>
//#include<networkit/graph/Graph.hpp>

bool temporalEdgeGreaterTimewise(const akt::TemporalEdge& lhs, const akt::TemporalEdge& rhs)
{
  return (lhs.when != rhs.when) ? (lhs.when > rhs.when)
    : ((lhs.from != rhs.from) ? (lhs.from < rhs.from)
       : (lhs.to < rhs.to));
}

// Stores the results of a bfs into a temporal graph: the shortest distance to it and the earliest possible arrival time there
struct VertexDistInfo
{
  // The shortest distance to the node
  int dist;
  // The earliest possible arrival time at the node
  int foremostTime;

};




// Struct for storing the helper arrays for each iteration of the outer loop of the algorithm for betweenness centrality
struct ShortestBetweennessData
{
  ShortestBetweennessData(int n, int T)
    : deltaDots{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
      sigmas{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
      preds{ std::vector<std::vector<std::unordered_set<VertexAppearance>>>(n, std::vector<std::unordered_set<VertexAppearance>>(T)) },
      dists{ std::vector<std::vector<int>>(n, std::vector<int>(T, -1)) },
      totalDists{ std::vector<int>(n, -1) },
      totalSigmas{ std::vector<int>(n, 0) },
      stack{ std::vector<VertexAppearance>() },
      foremostTimes{ std::vector<int>(n, -1) },
      deltaForemosts{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) }
  { }
  ShortestBetweennessData(const akt::Graph& g)
    : ShortestBetweennessData(g.N(), g.T())
  { }

  // The \delta_{s\cdot} array on each iteration
  std::vector<std::vector<double>> deltaDots;
  // Numbers of shortest paths from a source to each vertex appearance
  std::vector<std::vector<int>> sigmas;
  // Sets of predecessors for each vertex appearance
  std::vector<std::vector<std::unordered_set<VertexAppearance>>> preds;
  // Distances to each vertex appearance
  std::vector<std::vector<int>> dists;
  // Distances to each vertex (*not* vertex __appearance__)
  std::vector<int> totalDists;
  // Number of shortests paths from a source to each vertex (*not* vertex __appearance__)
  std::vector<int> totalSigmas;
  // Stack of vertex apperances in order of their discovery with the bfs
  std::vector<VertexAppearance> stack;


  // Additional arrays for handling shortest foremost computation
  // Earliest arrival times for each vertex
  std::vector<int> foremostTimes;
  // Delta array for foremost paths
  std::vector<std::vector<double>> deltaForemosts;
};



// Struct for storing the helper arrays for the prefix-foremost computation algorithm
struct PrefixBetweennessData
{
  PrefixBetweennessData(int n)
    : deltaDots(n, 1.0), sigmas(n, 0), preds(n), foremostTimes(n, -1)
  { }

  PrefixBetweennessData(const akt::Graph& g)
    : PrefixBetweennessData(g.N())
  { }

  // The \delta_{s\cdot} values for the foremost apperances
  std::vector<double> deltaDots;
  // Numbers of prefix-foremost paths from the source to each vertex
  std::vector<int> sigmas;
  // Sets of predecessors for each vertex
  std::vector<std::unordered_set<int>> preds;
  // foremost time to each vertex
  std::vector<int> foremostTimes;
  // Stack of nodes in order of their discovery with the "foremost-based" search
  std::vector<int> stack;
};
// (Re-) Initializes all the members in sbd for the next iteration of the outermost iteration of the shortest betwenness algorithm
void reinitializeHelperStructOptimal(const akt::Graph& g, int s, OptimalBetweennessData& sbd)
{
  // Reinitialize the appearance-based arrays (at the only relevant times)
  for (auto it = g.edges_cbegin(); it != g.edges_cend(); ++it) {
    auto te = *it;
    sbd.deltasvvt[te.to][te.when] = 0.0;
    sbd.sigmadot[te.to][te.when] = 0;
    sbd.pre[te.to][te.when].clear();
    sbd.totalSigmaT[te.to][te.when]=0;
    sbd.cur_best[te.to][te.when]=0;
    sbd.opt_walk[te.to][te.when]=Path(NULL,NULL);
  }
  for (int n = 0; n < g.N(); n++) {
    for (const auto& t: g.events) {
    sbd.deltadot[n][t] = 0.0;
    sbd.sigma[n][t] = 0;
    }
  }
  // Reinitialize the vertex-based arrays
  std::transform(sbd.optimalNode.cbegin(), sbd.optimalNode.cend(), sbd.optimalNode.begin(), [](auto i) { return std::numeric_limits<double>::infinity(); });
  std::transform(sbd.totalSigma.cbegin(), sbd.totalSigma.cend(), sbd.totalSigma.begin(), [](auto i) { return 0; });
  // Reinitialize the stack (though it should be empty already, so it's just a sanity operation)
  sbd.stack.clear();
  // Initialize the elements involving the source to appropriate values (if different from the defaults)



}

// (Re-) Initializes all the members in sbd for the next iteration of the outermost iteration of the shortest betwenness algorithm
void reinitializeHelperStruct(const akt::Graph& g, int s, ShortestBetweennessData& sbd)
{
  // Reinitialize the appearance-based arrays (at the only relevant times)
  for (auto it = g.edges_cbegin(); it != g.edges_cend(); ++it) {
    auto te = *it;
    sbd.deltaDots[te.to][te.when] = 0.0;
    sbd.deltaForemosts[te.to][te.when] = 0.0;
    sbd.sigmas[te.to][te.when] = 0;
    sbd.preds[te.to][te.when].clear();
    sbd.dists[te.to][te.when] = -1;
  }
  // Reinitialize the vertex-based arrays
  std::transform(sbd.totalDists.cbegin(), sbd.totalDists.cend(), sbd.totalDists.begin(), [](auto i) { return -1; });
  std::transform(sbd.totalSigmas.cbegin(), sbd.totalSigmas.cend(), sbd.totalSigmas.begin(), [](auto i) { return 0; });
  std::transform(sbd.foremostTimes.cbegin(), sbd.foremostTimes.cend(), sbd.foremostTimes.begin(), [](auto i) { return -1; });
  // Reinitialize the stack (though it should be empty already, so it's just a sanity operation)
  sbd.stack.clear();
  // Initialize the elements involving the source to appropriate values (if different from the defaults)
  sbd.sigmas[s][0] = 1;
  sbd.dists[s][0] = 0;
  sbd.totalDists[s] = 0;
  sbd.totalSigmas[s] = 1;
  sbd.foremostTimes[s] = 0;
}

// Computes all the distances to vertices, vertex appearances and the counts of corresponding shortest paths
void shortestComputeDistancesSigmas(const akt::Graph& g, bool strict, int s, ShortestBetweennessData& sbd)
{
  std::queue<VertexAppearance> q;
  q.push(VertexAppearance{ s, 0 });
  // BFS to find distances and sigma values
  while (!q.empty()) {
    auto cur = q.front();
    q.pop();
    // Go over all neighbours of the current vertex appearance
    for (int t = cur.time + strict; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep) {
      for (auto w : g.adjacencyList()[cur.v][t].neighbours) {
        if (sbd.dists[w][t] < 0) {
          sbd.dists[w][t] = sbd.dists[cur.v][cur.time] + 1;
          if (sbd.totalDists[w] < 0)
            sbd.totalDists[w] = sbd.dists[w][t];
          q.push(VertexAppearance{ w, t });
          sbd.stack.push_back(VertexAppearance{ w, t });
        }
        if (sbd.dists[w][t] == sbd.dists[cur.v][cur.time] + 1) {
          sbd.sigmas[w][t] += sbd.sigmas[cur.v][cur.time];
          sbd.preds[w][t].insert(cur);
          if (sbd.totalDists[w] == sbd.dists[cur.v][cur.time] + 1)
            sbd.totalSigmas[w] += sbd.sigmas[cur.v][cur.time];
        }
        if ((sbd.foremostTimes[w] < 0) || (t < sbd.foremostTimes[w])) {
          sbd.foremostTimes[w] = t;
        }
      }
    }
  }
}

void optimal_initialization(const akt::Graph& g, int s, OptimalBetweennessData& sbd, double (*cost)(Path, int, const akt::Graph&), FibonacciHeap<std::pair<double,std::pair<int,int>> > q, std::map<VertexAppearance, node<std::pair<double,std::pair<int,int>>>*> q_nodes)
{
  for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[s][t].nextTimestep) {
    if (g.adjacencyList()[s][t].neighbours.size() > 0)
      {
        sbd.opt_walk[s][t] = Path((std::vector<int>*) malloc(sizeof(std::vector<int>)), (std::vector<int>*) malloc(sizeof(std::vector<int>)));
        sbd.cur_best[s][t] = cost(sbd.opt_walk[s][t], t, g);
        //then nodes identifiers and times have to be positive
        sbd.pre[s][t].insert(VertexAppearance{ -1, -1 });
        std::pair<int,int> p_tmp (s,t);
        std::pair<double,std::pair<int,int>> p (sbd.cur_best[s][t], p_tmp);
        q_nodes[VertexAppearance{s,t}] = q.insert(p);
        sbd.optimalNode[s] = 0.0;

      }
  }
}

VertexAppearance pair_to_vertexappearance(std::pair<double, std::pair<int,int>> min_elem)
{
  auto pair_cur = min_elem.second;
  VertexAppearance cur;
  cur.v = pair_cur.first; cur.time = pair_cur.second;
  return cur;
}

void relax_resting(int b, int t, int tp, OptimalBetweennessData& sbd, FibonacciHeap<std::pair<double,std::pair<int,int>> > q, std::map<VertexAppearance, node<std::pair<double,std::pair<int,int>>>*> q_nodes, bool (*cmp)(double, double),double (*cost)(Path, int, const akt::Graph&), const akt::Graph& g)
{
  auto cnew = cost(sbd.opt_walk[b][t], tp, g);
  auto cold = cost(sbd.opt_walk[b][tp], tp, g);
  if (cmp(cnew, cold))
          {
            sbd.pre[b][tp].clear();
            sbd.cur_best[b][tp] = cnew;
            sbd.opt_walk[b][tp] = sbd.opt_walk[b][t];
            std::pair<int,int> p_tmp (b,tp);
            std::pair<double,std::pair<int,int>> p (cnew, p_tmp);
            if (q_nodes.count(VertexAppearance{b,tp}) == 1 )
              q.decreaseKey(q_nodes[VertexAppearance{b,tp}], p);
            else
              q_nodes[VertexAppearance{b,tp}] = q.insert(p);
            if (cnew < sbd.optimalNode[b])
              sbd.optimalNode[b] = cnew;

          }
}

void relax(int a, int b, int t, int tp, OptimalBetweennessData& sbd, FibonacciHeap<std::pair<double,std::pair<int,int>> > q,   std::map<VertexAppearance, node<std::pair<double,std::pair<int,int>>>*> q_nodes, bool (*cmp)(double, double),double (*cost)(Path, int, const akt::Graph&), const akt::Graph& g)
{
  if (sbd.pre[a][t].size() == 0)
    return;
  auto m = sbd.opt_walk[a][t];
  auto mp = m.clone();
  mp.add_link(a,b,tp);
  auto cnew = cost(mp, tp, g);
  auto cold = cost(sbd.opt_walk[b][tp], tp, g);
  if (cmp(cnew, cold))
          {
            sbd.pre[b][tp].clear();
            sbd.cur_best[b][tp] = cnew;
            sbd.opt_walk[b][tp] = mp;
            std::pair<int,int> p_tmp (b,tp);
            std::pair<double,std::pair<int,int>> p (cnew, p_tmp);
            if (q_nodes.count(VertexAppearance{b,tp}) == 1 )
              q.decreaseKey(q_nodes[VertexAppearance{b,tp}], p);
              else
                q_nodes[VertexAppearance{b,tp}] = q.insert(p);
            if (cnew < sbd.optimalNode[b])
              sbd.optimalNode[b] = cnew;

          }
  if (cnew == sbd.cur_best[b][tp])
    sbd.pre[b][tp].insert(VertexAppearance{a,t});
}

void optimalComputeDistancesSigmas(const akt::Graph& g, bool strict, int s, OptimalBetweennessData& sbd, double (*cost)(Path, int, const akt::Graph&), bool (*cmp)(double, double), std::string walk_type  )
{
  std::map<VertexAppearance, node<std::pair<double,std::pair<int,int>>>* > q_nodes;
  FibonacciHeap<std::pair<double,std::pair<int,int>> > q;
  optimal_initialization(g, s, sbd, cost, q, q_nodes);
  while (!q.isEmpty()) {
    auto min_elem = q.removeMinimum();
    VertexAppearance cur = pair_to_vertexappearance(min_elem);

    // Go over all neighbours of the current vertex appearance
    for (int t = cur.time + strict; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep_inv) {
      for (auto w : g.adjacencyList()[cur.v][t].neighbours_inv) {
        if (walk_type == "active")
          relax_resting(cur.v,cur.time,t, sbd, q, q_nodes, cmp, cost, g);
    }

    }

    for (int t = cur.time + strict; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep) {
      for (auto w : g.adjacencyList()[cur.v][t].neighbours_inv) {
        if (walk_type == "active")
          relax_resting(cur.v, cur.time, t, sbd, q, q_nodes, cmp, cost, g);
        relax(cur.v,w,cur.time,t,sbd,q,q_nodes,cmp,cost,g);
        }
      }
    }
}

// Empties the stack in sbd while updating the betweenness values of the nodes of the graph
void shortestEmptyStackUpdateBetweenness(int s, ShortestBetweennessData& sbd, std::vector<double>& betweenness, std::vector<double>& fmBetweenness)
{
  // Compute the delta array values (bottom-up) and update the betweenness values
  // Subtract the sum of the s-th row of the connectivity matrix to betweenness of s (follows from our formula)
  auto connectivityCorrection = std::accumulate(sbd.totalDists.cbegin(), sbd.totalDists.cend(), int(0),
                                                [](auto acc, auto dist) { return acc + (dist >= 0); });
  betweenness[s] -= connectivityCorrection;
  fmBetweenness[s] -= connectivityCorrection;
  while (!sbd.stack.empty()) {
    auto cur = sbd.stack.back();
    sbd.stack.pop_back();
    // Handle "base-cases," i.e. times where the appearance is the end of a shortest (foremost) path to the vertex in question
    if (sbd.dists[cur.v][cur.time] == sbd.totalDists[cur.v])
      sbd.deltaDots[cur.v][cur.time] += static_cast<double>(sbd.sigmas[cur.v][cur.time])
        / static_cast<double>(sbd.totalSigmas[cur.v]);
    if (cur.time == sbd.foremostTimes[cur.v])
      sbd.deltaForemosts[cur.v][cur.time] += 1.0;
    for (auto pred : sbd.preds[cur.v][cur.time]) {
      // Compute the currently considered summand of the sum over successors of pred and update betweenness centrality and deltaDots
      const auto shortestSummand = (static_cast<double>(sbd.sigmas[pred.v][pred.time])
                                    / static_cast<double>(sbd.sigmas[cur.v][cur.time]))
        * sbd.deltaDots[cur.v][cur.time];
      sbd.deltaDots[pred.v][pred.time] += shortestSummand;
      betweenness[pred.v] += shortestSummand;
      const auto foremostSummand = (static_cast<double>(sbd.sigmas[pred.v][pred.time])
                                    / static_cast<double>(sbd.sigmas[cur.v][cur.time]))
        * sbd.deltaForemosts[cur.v][cur.time];
      sbd.deltaForemosts[pred.v][pred.time] += foremostSummand;
      fmBetweenness[pred.v] += foremostSummand;
    }
  }
}

void UpdateBetweenness(OptimalBetweennessData& sbd, int n, const std::vector<int> &events)
{
  for (int i = 0; i < n; i++)
    {
      for (const auto& t: events) {
        sbd.betweenness[i][t] = sbd.betweenness[i][t] + sbd.deltadot[i][t];
      }
    }
}

void optimalUpdateBetweenness(int s, const akt::Graph& g, OptimalBetweennessData& sbd, double (*cost)(Path, int, const akt::Graph&), bool (*cmp)(double, double), std::string walk_type)
{
  std::vector<std::vector<double>> betweenness;
  Predecessor G = PredecessorGraph(g, sbd.pre, s);
  std::pair<std::unordered_set<int>, std::unordered_set<int>> p;
  p = RemoveInfiniteFromPredecessor(s, G, sbd, cost, cmp, walk_type, g);
  VolumePathAt(G, s, sbd);
  OptimalSigma(s, G, sbd, g, cost);
  ComputeDeltaSvvt(G, s, sbd);
  PredecessorGraphToOrdered(G, g.maximalTimestep(), g.events_rev);
  std::map<VertexAppearance, VertexAppearance> preced;
  GeneralContribution(g, G, s, sbd, preced, walk_type);
  UpdateBetweenness(sbd, g.N(),  g.events);
}

void prefixDoForemostBasedSearch(const akt::Graph& g, int s, PrefixBetweennessData& pbd)
{
  auto q = std::priority_queue<akt::TemporalEdge, std::vector<akt::TemporalEdge>, decltype(&temporalEdgeGreaterTimewise)>(&temporalEdgeGreaterTimewise);
  for (auto it = g.firstEdgeFromVertex(s); (it != g.edges_cend()) && ((*it).from == s); ++it)
    q.push(*it);
  while (!q.empty()) {
    auto te = q.top();
    q.pop();
    if (pbd.foremostTimes[te.to] < 0) {
      pbd.foremostTimes[te.to] = te.when;
      pbd.stack.push_back(te.to);
      if (te.when < g.maximalTimestep())
        for (auto it = g.firstEdgeFromAppearance(te.to, te.when + 1); (it != g.edges_cend()) && ((*it).from == te.to); ++it)
          q.push(*it);
    }
    if (te.when == pbd.foremostTimes[te.to]) {
      pbd.sigmas[te.to] += pbd.sigmas[te.from];
      pbd.preds[te.to].insert(te.from);
    }
  }
}

void prefixEmptyStackUpdateBetweenness(int s, PrefixBetweennessData& pbd, std::vector<double>& betweenness)
{
  betweenness[s] -= std::accumulate(pbd.foremostTimes.cbegin(), pbd.foremostTimes.cend(), int(0),
                                    [](auto acc, auto time) { return acc + (time >= 0); });
  while (!(pbd.stack.empty())) {
    auto w = pbd.stack.back();
    pbd.stack.pop_back();
    for (auto v : pbd.preds[w]) {
      // Compute the currently considered summand of the sum over successors of pred and update betweenness centrality and deltaDots
      auto summand = (static_cast<double>(pbd.sigmas[v])
                      / static_cast<double>(pbd.sigmas[w]))
        * pbd.deltaDots[w];
      pbd.deltaDots[v] += summand;
      betweenness[v] += summand;
    }
  }
}

namespace akt {
  bool smaller(double x, double y)
  {
    if (x <= y)
      return true;
    return false;
  }
  // Computes the betweenness measures
  std::vector<std::vector<double>> optimalBetweenness(const Graph& g, bool strict, std::string cost, std::string cmp, std::string walk_type)
  {
    double (*cost2)(Path, int, const akt::Graph&);
    bool (*cmp2)(double, double);

    if (cost == "shortestfastest")
      cost2 = &co_sfp;
    if (cmp == "le")
      cmp2 = &smaller;
    auto sbd = OptimalBetweennessData(g);
    // Stores the betweenness values; initialize to 1 because of the formula for betweenness having a constant +1
    for (int s = 0; s < g.N(); ++s) {
      reinitializeHelperStructOptimal(g, s, sbd);
      optimalComputeDistancesSigmas(g, strict, s, sbd, cost2, cmp2, walk_type);
      optimalUpdateBetweenness(s, g, sbd, cost2, cmp2, walk_type);
    }
    return sbd.betweenness;
  }

  // Computes the betweenness measures
  std::pair<std::vector<double>, std::vector<double>> shortestBetweenness(const Graph& g, bool strict)
  {
    auto sbd = ShortestBetweennessData(g);
    // Stores the betweenness values; initialize to 1 because of the formula for betweenness having a constant +1
    auto shortestBetweenness = std::vector<double>(g.N(), 1.0);
    auto foremostBetweenness = std::vector<double>(g.N(), 1.0);
    for (int s = 0; s < g.N(); ++s) {
      reinitializeHelperStruct(g, s, sbd);
      shortestComputeDistancesSigmas(g, strict, s, sbd);
      shortestEmptyStackUpdateBetweenness(s, sbd, shortestBetweenness, foremostBetweenness);
    }
    return { shortestBetweenness, foremostBetweenness };
  }

  std::vector<double> prefixForemostBetweenness(const Graph& g)
  {
    auto res = std::vector<double>(g.N(), 1.0);
    for (int s = 0; s < g.N(); ++s) {
      auto pbd = PrefixBetweennessData(g);
      pbd.foremostTimes[s] = 0;
      pbd.sigmas[s] = 1;
      prefixDoForemostBasedSearch(g, s, pbd);
      prefixEmptyStackUpdateBetweenness(s, pbd, res);
    }

    return res;
  }


  // Returns the shortest betweenness for a static graph (i.e. all edges have timestamp 0) using Brandes' algorithm
  std::vector<double> shortestBetweennessStatic(const Graph& g)
  {
    auto res = std::vector<double>(g.N(), 0);
    for (int s = 0; s < g.N(); ++s) {
      auto stack = std::vector<int>();
      auto preds = std::vector<std::vector<int>>(g.N(), std::vector<int>());
      auto sigma = std::vector<int>(g.N(), 0);
      auto dists = std::vector<int>(g.N(), -1);
      sigma[s] = 1;
      dists[s] = 0;
      std::queue<int> q;
      q.push(s);
      while (!q.empty()) {
        auto v = q.front();
        q.pop();
        stack.push_back(v);
        for (auto w : (g.adjacencyList()[v][0]).neighbours) {
          if (dists[w] < 0) {
            dists[w] = dists[v] + 1;
            q.push(w);
          }
          if (dists[w] == (dists[v] + 1)) {
            sigma[w] += sigma[v];
            preds[w].push_back(v);
          }
        }
      }
      auto delta = std::vector<double>(g.N(), 0.0);
      while (!stack.empty()) {
        auto w = stack.back();
        stack.pop_back();
        for (auto v : preds[w])
          delta[v] += static_cast<double>(sigma[v]) / static_cast<double>(sigma[w])* (1 + delta[w]);
        if (w != s)
          res[w] += delta[w];
      }
    }

    return res;
  }

} // end namespace akt
