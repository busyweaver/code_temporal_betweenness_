#include <boost/heap/fibonacci_heap.hpp>
#include "predecessor_graph.h"
//#include "paths.h"
#include <numeric>
#include <queue>
//#include<networkit/graph/Graph.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include<random>
#include<iostream>
using namespace std;
using namespace boost::heap;


struct node
{
  pair<double,pair<int,int>> id;

  node(pair<double,pair<int,int>> i)
    : id(i)
  { }
};

struct compare_node
{
  bool operator()(const node& n1, const node& n2) const
  {
    auto a = n1.id;
    auto b = n2.id;
    if (a.first > b.first)
      {
        return true;
      }
    else if (a.first == b.first)
      {
        if (a.second.first > b.second.first)
          return true;
        else if(a.second.first == b.second.first)
          return (a.second.second > b.second.second);
        else
          return false;
      }
    else
      return false;
  }
};


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
void reinitializeHelperStructOptimal(const akt::Graph& g, int s, OptimalBetweennessData& sbd, std::vector<long int> &vis)
{
  // Reinitialize the appearance-based arrays (at the only relevant times)
  // printf("okok\n");
  int T = g.events.size();
  while (!vis.empty()) {
    auto el = vis.back();
    vis.pop_back();
    long int n = el / T;
    long int t = el % T;
    //    std::cout << "init " << n << " " << t << "\n";
    // std::cout <<" T ->" << T;
    // for (int n = 0; n < g.N(); n++) {
    //   for (int t = 0; t<T ; t++) {
    sbd.deltasvvt[n][t] = 0.0;
    sbd.deltadot[n][t] = 0.0;
    if(sbd.opt_walk[n][t] != nullptr)
      {
        if (sbd.opt_walk[n][t]->arrival() == t)
          delete sbd.opt_walk[n][t];
      }

    sbd.opt_walk[n][t]=nullptr;
    sbd.pre[n][t].clear();
    sbd.cur_best[n][t]=std::numeric_limits<double>::infinity();

    sbd.totalSigmaT[n][t]=0;
    sbd.sigma[n][t] = 0;
    sbd.sigmadot[n][t] = 0;
    delete sbd.opt_walk[n][t];
  }

  for (int n = 0; n < g.N(); n++) {
    sbd.totalSigma[n] = 0;
    sbd.optimalNode[n] = std::numeric_limits<double>::infinity();
  }

}

void reinitializeHelperStructOptimalBoost(const akt::Graph& g, int s, OptimalBetweennessData& sbd)
{
  // Reinitialize the appearance-based arrays (at the only relevant times)
  int T = g.events.size();
  while (!sbd.visited.empty()) {
    auto el = sbd.visited.back();
    sbd.visited.pop_back();
    long int n = el / T;
    long int t = el % T;
    sbd.deltasvvt[n][t] = 0.0;
    sbd.deltadot[n][t] = 0.0;
    sbd.pre[n][t].clear();
    sbd.cur_best[n][t]=std::numeric_limits<double>::infinity();
    sbd.totalSigmaT[n][t]=0;
    sbd.sigma[n][t] = 0;
    sbd.sigmadot[n][t] = 0;
    delete sbd.opt_walk[n][t];
  }
  for (int n = 0; n < g.N(); n++) {
    sbd.totalSigma[n] = 0;
    sbd.optimalNode[n] = std::numeric_limits<double>::infinity();
  }
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


void display_tot(OptimalBetweennessData& sbd)
{


  printf("pre\n");
  for(auto &elem : sbd.pre)
    {
      for(auto &elem2 : elem.second)
        {
          std::cout<<"pre  " <<  elem.first << " " << elem2.first << "\n";
          for(auto &elem3 : elem2.second)
            {
              std::cout<<"::    " <<  elem3.v << " " << elem3.time;
            }
          std::cout << "\n";
        }
    }

  printf("cur_best\n");
  for (int i = 0; i < sbd.cur_best.size(); i++)
    {
      for (int j = 0; j < sbd.cur_best[i].size(); j++)
        {
          std::cout<<" " <<sbd.cur_best[i].at(j);
        }
      std::cout<< " " << "\n";
    }


  printf("sigmadot\n");
for (int i = 0; i < sbd.cur_best.size(); i++)
  {
    for (int j = 0; j < sbd.cur_best[i].size(); j++)
      {
        std::cout<<" " <<sbd.sigmadot[i].at(j); 
      }
    std::cout<< " " << "\n";
  }
 printf("sigmatotT\n");
 for (int i = 0; i < sbd.totalSigmaT.size(); i++)
   {
     for (int j = 0; j < sbd.totalSigmaT[i].size(); j++)
       {
         std::cout<<" " <<sbd.totalSigmaT[i].at(j); 
       }
     std::cout<< " " << "\n";
   }

 printf("sigma\n");
 for (int i = 0; i < sbd.sigma.size(); i++)
   {
     for (int j = 0; j < sbd.sigma[i].size(); j++)
       {
         std::cout<<" " <<sbd.sigma[i].at(j); 
       }
     std::cout<< " " << "\n";
   }

 printf("deltasvvt\n");
 for (int i = 0; i < sbd.deltasvvt.size(); i++)
   {
     for (int j = 0; j < sbd.deltasvvt[i].size(); j++)
       {
         std::cout<<" " <<sbd.deltasvvt[i].at(j); 
       }
     std::cout<< " " << "\n";
   }

 printf("deltadot\n");
 for (int i = 0; i < sbd.deltadot.size(); i++)
   {
     for (int j = 0; j < sbd.deltadot[i].size(); j++)
       {
         std::cout<<" " <<sbd.deltadot[i].at(j); 
       }
     std::cout<< " " << "\n";
   }
 printf("betweenness\n");
 for (int i = 0; i < sbd.betweenness.size(); i++)
   {
     for (int j = 0; j < sbd.betweenness[i].size(); j++)
       {
         std::cout<<" " <<sbd.betweenness[i].at(j); 
       }
     std::cout<< " " << "\n";
   }
 printf("betweenness_exact\n");
 for (int i = 0; i < sbd.betweenness_exact.size(); i++)
   {
     for (int j = 0; j < sbd.betweenness_exact[i].size(); j++)
       {
         std::cout<<" " <<sbd.betweenness_exact[i].at(j); 
       }
     std::cout<< " " << "\n";
   }

 printf("optimal_node\n");
 for (int i = 0; i < sbd.optimalNode.size(); i++)
   {
     std::cout<<" " <<sbd.optimalNode[i];
   }


 printf("\nsigma_total\n");
 for (int i = 0; i < sbd.totalSigma.size(); i++)
   {
     std::cout<<" " <<sbd.totalSigma[i];
   }


 printf("\nbetweenness_exact_sum\n");
 for (int i = 0; i < sbd.betweenness_exact.size(); i++)
   {
     double s = 0;
     for (int j = 0; j < sbd.betweenness_exact[i].size(); j++)
       {
         s = s + sbd.betweenness_exact[i].at(j);
       }
     std::cout<< s <<" " << "\n";
   }
 printf("end display \n");

}



void dijkstra_initialization(const akt::Graph& g, int s, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph&), fibonacci_heap<node, boost::heap::compare<compare_node>> &q, std::vector<long int> &visited)
{
  int T = g.events.size();
  for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[s][t].nextTimestep) {

    if (g.adjacencyList()[s][t].neighbours.size() > 0)
      {
      Path* pp = new Path(nullptr, s, g.events[t]);
      sbd.opt_walk[s][t] = pp;
      sbd.cur_best[s][t] = cost(sbd.opt_walk[s][t], g.events[t], g);
      sbd.pre[s][t].insert(VertexAppearance{ s, -1 });
      //add to visited since not added to heap, update : not sur the following line is needed
      visited.push_back(s*T + t);
      q.push(node({ 0.0 , {s,t}  }));
      }
    }
  sbd.optimalNode[s] = 0.0;

}

void optimalInitializationBoost(const akt::Graph& g, int s, OptimalBetweennessData& sbd, std::queue<int> &q, double untilTime)
{
  int T = g.events.size();
  for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()) && (g.events[t] <= untilTime); t = g.adjacencyList()[s][t].nextTimestep) {
    if (g.adjacencyList()[s][t].neighbours.size() > 0)
      {
        //        Path* pp = new Path(nullptr, s, g.events[t]);
        sbd.opt_walk[s][t] = nullptr;
        sbd.cur_best[s][t] = 0;
        sbd.pre[s][t].insert(VertexAppearance{ s, -1 });
        sbd.optimalNode[s] = 0.0;
        sbd.sigmadot[s][t] = 1;
        q.push(s*T + t);
      }
    sbd.totalSigma[s] = 1;
    }
}

VertexAppearance pair_to_vertexappearance(std::pair<double, std::pair<int,int>> min_elem)
{
  auto pair_cur = min_elem.second;
  VertexAppearance cur;
  cur.v = pair_cur.first; cur.time = pair_cur.second;
  return cur;
}

void relax_resting(int b, int t, int tp, OptimalBetweennessData& sbd, fibonacci_heap<node, compare<compare_node>> &q, map<pair<double,pair<int,int>>,fibonacci_heap<node, compare<compare_node>>::handle_type> ma,bool (*cmp)(double, double),double (*cost)(Path*, int, const akt::Graph&), const akt::Graph& g)
{
  int T = g.events.size();
  auto cnew = cost(sbd.opt_walk[b][t], g.events[tp], g);
  auto cold = cost(sbd.opt_walk[b][tp], g.events[tp], g);
  if (cmp(cnew, cold))
          {
            sbd.pre[b][tp].clear();
            sbd.cur_best[b][tp] = cnew;
            sbd.opt_walk[b][tp] = sbd.opt_walk[b][t];
            std::pair<int,int> p_tmp;
            p_tmp.first = b;
            p_tmp.second = tp;
            std::pair<double,std::pair<int,int>> p;
            p.first = cnew;
            p.second = p_tmp;
            std::pair<int,int> p_tmp2;
            p_tmp2.first = b;
            p_tmp2.second = tp;
            std::pair<double,std::pair<int,int>> p2;
            p2.first = cold;
            p2.second = p_tmp2;
            if (ma.count(p2) == 1 )
              {
                q.increase(ma[p2], p);
                //                std::cout << "DECREASE.................. \n"; 
              }
            else
              q.push(node(p));
            if (cnew < sbd.optimalNode[b])
              sbd.optimalNode[b] = cnew;
          }
}

void relax(int a, int b, int t, int tp, OptimalBetweennessData& sbd, fibonacci_heap<node, boost::heap::compare<compare_node>> &q, map<pair<double,std::pair<int,int>>,fibonacci_heap<node, compare<compare_node>>::handle_type> ma, bool (*cmp)(double, double),double (*cost)(Path*, int, const akt::Graph&), const akt::Graph& g)
{
  int T = g.events.size();
  if (sbd.pre[a][t].size() == 0)
    return;
  auto m = sbd.opt_walk[a][t];
  auto mp = new Path(m,b,g.events[tp]);
  auto cnew = cost(mp, g.events[tp], g);
  auto cold = cost(sbd.opt_walk[b][tp], g.events[tp], g);
  if (cmp(cnew, cold)){
    auto tmp = sbd.opt_walk[b][tp];
            sbd.pre[b][tp].clear();
            sbd.cur_best[b][tp] = cnew;
            sbd.opt_walk[b][tp] = mp;
            std::pair<int,int> p_tmp;
            p_tmp.first = b;
            p_tmp.second = tp;
            std::pair<double,std::pair<int,int>> p;
            p.first = cnew;
            p.second = p_tmp;
            std::pair<int,int> p_tmp2;
            p_tmp2.first = b;
            p_tmp2.second = tp;
            std::pair<double,std::pair<int,int>> p2;
            p2.first = cold;
            p2.second = p_tmp2;
            //printf("optimal\n");
            if (ma.count(p2) == 1 ){
              // in reality we decrease 
              q.increase(ma[p2], p);
              //display(tmp);
              //display(mp);
              //std::cout <<"old " << cold << " " << p.first << " " << p.second.first << " " << p.second.second  << "DECREASE.................. \n"; 
              }
            else{
              q.push(node(p));
              }
            if (cnew < sbd.optimalNode[b])
              {
                sbd.optimalNode[b] = cnew;
              }

          }
  else
    delete mp;
  if (cnew == sbd.cur_best[b][tp] and cnew != std::numeric_limits<double>::infinity())
    {
      sbd.pre[b][tp].insert(VertexAppearance{a,t});
      // sbd.sigmadot[b][tp] += sbd.sigmadot[a][t];
      // if(cnew == sbd.optimalNode[b])
      //   sbd.totalSigma[b] = sbd.totalSigma[b] + sbd.sigmadot[a][t];
    }
  //printf("fin\n");

}

/* ************************* bellman-ford non-optimal **************************  */
void relax_resting_bf(int b, int t, int tp, OptimalBetweennessData& sbd, bool (*cmp)(double, double),double (*cost)(Path*, int, const akt::Graph&), const akt::Graph& g, std::vector<long int>& vis)
{
  int T = g.events.size();
  auto cnew = cost(sbd.opt_walk[b][t], g.events[tp], g);
  auto cold = cost(sbd.opt_walk[b][tp], g.events[tp], g);
  if (cmp(cnew, cold)){
    vis.push_back(b*T + tp);
    sbd.pre[b][tp].clear();
    sbd.cur_best[b][tp] = cnew;
    sbd.opt_walk[b][tp] = sbd.opt_walk[b][t];
    
  }
}

void relax_bf(int a, int b, int t, int tp, OptimalBetweennessData& sbd, bool (*cmp)(double, double),double (*cost)(Path*, int, const akt::Graph&), const akt::Graph& g, std::vector<long int>& vis, std::string walk_type, int s)
{
  if (sbd.pre[a][t].size() == 0)
    return;
  int T = g.events.size();
  auto m = sbd.opt_walk[a][t];
  auto mp = new Path(m,b,g.events[tp]);
  auto cnew = cost(mp, g.events[tp], g);
  auto cold = cost(sbd.opt_walk[b][tp], g.events[tp], g);
  if (cmp(cnew, cold)){
    vis.push_back(b*T + tp);
    auto tmp = sbd.opt_walk[b][tp];
    sbd.pre[b][tp].clear();
    sbd.cur_best[b][tp] = cnew;
    sbd.opt_walk[b][tp] = mp;
    std::pair<int,int> p_tmp;
    p_tmp.first = b;
    p_tmp.second = tp;
    std::pair<double,std::pair<int,int>> p;
    p.first = cnew;
    p.second = p_tmp;
    std::pair<int,int> p_tmp2;
    p_tmp2.first = b;
    p_tmp2.second = tp;
    std::pair<double,std::pair<int,int>> p2;
    p2.first = cold;
    p2.second = p_tmp2;
    if (cnew < sbd.optimalNode[b])
      sbd.optimalNode[b] = cnew;
    if (walk_type == "active")
      {
        for (int tpp = tp+1; (tpp >= 0) && (tpp <= g.maximalTimestep()); tpp = g.adjacencyList()[b][tpp].nextTimestep_inv) {
          for (auto c : g.adjacencyList()[b][tpp].neighbours_inv) {
            if (b != s)
              {
                relax_resting_bf(b,tp,tpp, sbd, cmp, cost, g, vis);
              }

          }
        }
      }
  }
  else
    delete mp;
  if (cnew == sbd.cur_best[b][tp] and cnew != std::numeric_limits<double>::infinity())
    {

      sbd.pre[b][tp].insert(VertexAppearance{a,t});
      // sbd.sigmadot[b][tp] += sbd.sigmadot[a][t];
      // if(cnew == sbd.optimalNode[b])
      //   sbd.totalSigma[b] = sbd.totalSigma[b] + sbd.sigmadot[a][t];
    }

}

//Bellman-Ford variant see article Rymar et al. Towards Classifying the Polynomial-Time Solvability of Temporal Betweenness Centrality
std::vector<long int>  BF(const akt::Graph& g, int s, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph&), bool (*cmp)(double, double), std::string walk_type)
{
  int T = g.events.size();
  std::vector<long int> visited;
  //initialization step
  for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[s][t].nextTimestep) {

    if (g.adjacencyList()[s][t].neighbours.size() > 0)
      {
        Path* pp = new Path(nullptr, s, g.events[t]);
        sbd.opt_walk[s][t] = pp;
        sbd.cur_best[s][t] = cost(sbd.opt_walk[s][t], g.events[t], g);
        sbd.pre[s][t].insert(VertexAppearance{ s, -1 });
        visited.push_back(s*T + t);
      }
  }
  sbd.optimalNode[s] = 0.0;
  //rest of the algorithm
  for(int i = 0; i < g.N()*g.T(); i++)
    {
      for(int v = 0; v < g.N(); v++)
        {
          for (int tp = g.minimalTimestep(); (tp >= 0) && (tp <= g.maximalTimestep()); tp = g.adjacencyList()[v][tp].nextTimestep)
            for (auto w : g.adjacencyList()[v][tp].neighbours) {
              for(int t = g.minimalTimestep(); t <= tp; t++)
                if((! (v == s && t != tp)) ){
                  relax_bf(v,w,t,tp,sbd,cmp,cost,g,visited, walk_type, s);
                }
            }
        }
    }

  return visited;
}
/* ************************* end bellman-ford non-optimal **************************  */



/* ************************* start BFS  **************************  */
std::vector<long int> BFS(const akt::Graph& g, int s, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph&), bool (*cmp)(double, double), std::string& walk_type)
{
  int T = g.events.size();
  std::queue<int>* q = new queue<int>();
  std::queue<int>* qp = new queue<int>();
  std::vector<long int> visited;
  // int l = 0;

  // initialization step
  for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[s][t].nextTimestep) {
    if (g.adjacencyList()[s][t].neighbours.size() > 0)
      {
        //        Path* pp = new Path(nullptr, s, g.events[t]);
        sbd.opt_walk[s][t] = nullptr;
        sbd.cur_best[s][t] = 0.0;
        sbd.pre[s][t].insert(VertexAppearance{ s, -1 });
        visited.push_back(s*T + t);

        //sbd.sigmadot[s][t] = 1;
        (*q).push(s*T + t);
      }
    // sbd.totalSigma[s] = 1;
  }
  sbd.optimalNode[s] = 0.0;
  // end initialization step
  while ((*q).size() > 0 ) {
    while((*q).size() > 0){
      auto elem = (*q).front();
      (*q).pop();
      auto curv = elem / T;
      auto curtime = elem % T;
      visited.push_back(curv*T + curtime);
      for (int t = curtime; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[curv][t].nextTimestep) {
        for (auto w : g.adjacencyList()[curv][t].neighbours) {
          if((! (curv == s && t != curtime)) ){
            auto m = sbd.opt_walk[curv][curtime];
            auto mp = new Path(m,w,g.events[t]);
            auto cnew = cost(mp, g.events[t], g);
            if (sbd.cur_best[w][t] == std::numeric_limits<double>::infinity() || sbd.pre[w][t].size() == 0 ){
              sbd.cur_best[w][t] = cnew;
              sbd.opt_walk[w][t] = mp;
              if (cnew < sbd.optimalNode[w])
                sbd.optimalNode[w] = cnew;
              (*qp).push(w*T + t);
              if (walk_type == "active")
                {
                  for (int tpp = t; (tpp >= 0) && (tpp <= g.maximalTimestep()); tpp = g.adjacencyList()[w][tpp].nextTimestep_inv) {
                    auto cnew = cost(sbd.opt_walk[w][t], g.events[tpp], g);
                    if (sbd.cur_best[w][tpp] > cnew){
                      sbd.visited.push_back(w*T + tpp);
                      //sbd.pre[w][tpp].clear();
                      sbd.cur_best[w][tpp] = cnew;
                    }
                  }
                }
            }
                  else
                    delete mp;

            if (sbd.cur_best[w][t] == cnew){
              sbd.pre[w][t].insert(VertexAppearance{curv,curtime});
              // sbd.sigmadot[w][t] += sbd.sigmadot[curv][curtime];
              // if(cnew == sbd.optimalNode[w])
              //   sbd.totalSigma[w] = sbd.totalSigma[w] + sbd.sigmadot[curv][curtime];
            }
          }
        }
      }
    }
    auto tmp = q;
    q = qp;
    qp = tmp;
    //l = l + 1;
  }
  delete q;
  delete qp;
  return visited;
}



/* ************************* end BFS  **************************  */

std::vector<long int>  dijkstra(const akt::Graph& g, int s, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph&), bool (*cmp)(double, double), std::string walk_type)
{

  int T = g.events.size();
  map<pair<double,std::pair<int,int>>,fibonacci_heap<node, compare<compare_node>>::handle_type> ma;
  fibonacci_heap<node, boost::heap::compare<compare_node>> q;
  std::vector<long int> visited;
  dijkstra_initialization(g, s, sbd, cost, q, visited);
  //printf("normalement boucle après init %ld,\n",q.size());
  int j = 0;
  while (!q.empty()) {

    auto min_elem = q.top();
    q.pop();
    ma.erase(min_elem.id);
    VertexAppearance cur = pair_to_vertexappearance(min_elem.id);
    visited.push_back(cur.v*T + cur.time);
    if(cur.v != s && walk_type == "active")
      {
        for (int t = cur.time; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep_inv) {
          for (auto w : g.adjacencyList()[cur.v][t].neighbours_inv) {
            relax_resting(cur.v,cur.time,t, sbd, q, ma, cmp, cost, g);
          }
        }
      }

    for (int t = cur.time; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep) {
      //std::cout << "nei time " << t << "cur.v " << cur.v << "cur.t " << cur.time <<  "\n"    ;
      for (auto w : g.adjacencyList()[cur.v][t].neighbours) {
        //std::cout << "nei time " << t << "w" << w << "\n";
        if((! (cur.v == s && t != cur.time)) ){
          // next condition could be removed
          if((cur.v == s || t >= (cur.time)) ){
            //std::cout << "nei time " << t << "cur.v " << cur.v << "cur.t " << cur.time << "w "<< w<<  "\n"    ;
            relax(cur.v,w,cur.time,t,sbd,q,ma,cmp,cost,g);
          }
        }
      }
    }
  }

  return visited;
}

void optimalComputeDistancesSigmasBoost(const akt::Graph& g, bool stri, int s, OptimalBetweennessData& sbd, std::string& walk_type, double k, int foremost, double untilTime)
{
  int T = g.events.size();
  int strict = 0;
  if(stri)
    strict = 1;
  std::queue<int>* q = new queue<int>();
  std::queue<int>* qp = new queue<int>();
  int l = 0;
  optimalInitializationBoost(g, s, sbd, *q, untilTime);
  while ((*q).size() > 0 ) {
    while((*q).size() > 0){
      auto elem = (*q).front();
      (*q).pop();
      auto curv = elem / T;
      auto curtime = elem % T;
      sbd.visited.push_back(curv*T + curtime);
      for (int t = curtime; (t >= 0) && (t <= g.maximalTimestep()) && (g.events[t] <= untilTime); t = g.adjacencyList()[curv][t].nextTimestep) {
        for (auto w : g.adjacencyList()[curv][t].neighbours) {
          if((! (curv == s && t != curtime)) ){
            if((curv == s || t >= (curtime+strict)) && (g.events[t] - g.events[(curtime+strict)] <= k ) ){
                  if (sbd.cur_best[w][t] == std::numeric_limits<double>::infinity() || (sbd.cur_best[w][t] >= l+1 && sbd.pre[w][t].size() == 0) ){
                    sbd.cur_best[w][t] = l+1;
                    if (!foremost)
                      {
                        if (l+1 < sbd.optimalNode[w])
                          sbd.optimalNode[w] = l+1; 
                      }
                    else
                      {
                        if (t < sbd.optimalNode[w])
                          {
                            sbd.optimalNode[w] = t;
                            sbd.totalSigma[w] = 0;
                          }

                      }
                    (*qp).push(w*T + t);
                    if (walk_type == "active")
                      {
                        for (int tpp = t + strict; (tpp >= 0) && (tpp <= g.maximalTimestep()) && (g.events[tpp] - g.events[(t+strict)] <= k ) && (g.events[tpp] <= untilTime); tpp = g.adjacencyList()[w][tpp].nextTimestep_inv) {
                          if (sbd.cur_best[w][tpp] > l+1){
                            sbd.visited.push_back(w*T + tpp);
                            //sbd.pre[w][tpp].clear();
                            sbd.cur_best[w][tpp] = l+1;
                          }
                        }
                      }
                  }
                  if (sbd.cur_best[w][t] == l+1){
                      sbd.pre[w][t].insert(VertexAppearance{curv,curtime});
                      sbd.sigmadot[w][t] = sbd.sigmadot[w][t] + sbd.sigmadot[curv][curtime];
                      if(!foremost)
                        {
                          if (l+1 == sbd.optimalNode[w])
                            sbd.totalSigma[w] = sbd.totalSigma[w] + sbd.sigmadot[curv][curtime];
                        }
                      else
                        {
                          if (t == sbd.optimalNode[w])
                            {
                              sbd.totalSigma[w] = sbd.totalSigma[w] + sbd.sigmadot[curv][curtime];
                            }

                        }
                  }
              }
          }
        }
      }
    }
    auto tmp = q;
    q = qp;
    qp = tmp;
    l = l + 1;
  }
  delete q;
  delete qp;
}

// main algorithm for shortest paths variants on passive case
void shortestEmptyStackUpdateBetweennessBoost(int s, OptimalBetweennessData& sbd, const akt::Graph& g, int foremost)
{
  int T = g.events.size();
  for (auto it = sbd.visited.rbegin(); it != sbd.visited.rend(); ++it){
    //  while (!sbd.visited.empty()) {
    auto tmp = *it;
    auto curv = tmp / T;
    auto curtime = tmp % T;
    //sbd.visited.pop_back();
    // Handle "base-cases," i.e. times where the appearance is the end of a shortest (foremost) path to the vertex in question
    if(!foremost)
      {
        if (sbd.cur_best[curv][curtime] == sbd.optimalNode[curv])
          sbd.deltadot[curv][curtime] += static_cast<double>(sbd.sigmadot[curv][curtime])
            / static_cast<double>(sbd.totalSigma[curv]);
      }
    else
      {
        if (curtime == sbd.optimalNode[curv])
          sbd.deltadot[curv][curtime] += static_cast<double>(sbd.sigmadot[curv][curtime])
            / static_cast<double>(sbd.totalSigma[curv]);
      }

    for (auto pred : sbd.pre[curv][curtime]) {
      if(pred.time != -1){
      // Compute the currently considered summand of the sum over successors of pred and update betweenness centrality and deltaDots
      const auto shortestSummand = (static_cast<double>(sbd.sigmadot[pred.v][pred.time])
                                    / static_cast<double>(sbd.sigmadot[curv][curtime]))
        * sbd.deltadot[curv][curtime];
      sbd.deltadot[pred.v][pred.time] += shortestSummand;
      if(pred.v != s)
        sbd.betweenness_exact[pred.v][pred.time] += shortestSummand;
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

void UpdateBetweenness(OptimalBetweennessData& sbd, int T, std::unordered_set<long int>& visited)
{
  for(auto &el : visited)
    {
      auto i = el/T;
      auto t = el%T;
      sbd.betweenness[i][t] = sbd.betweenness[i][t] + sbd.deltadot[i][t];
      sbd.betweenness_exact[i][t] = sbd.betweenness_exact[i][t] + sbd.deltadot[i][t];
    }
}

void totalBetweenness_compute(OptimalBetweennessData& sbd, int n, int T)
{
  for (int i = 0; i < n; i++)
    {
      double s = 0;
      for (int t = 0; t<T;t++)
        s += sbd.betweenness_exact[i][t];
      sbd.totalBetweenness[i] = s;
    }
}

void UpdateBetweenness_exact(OptimalBetweennessData& sbd, int T, int s, std::unordered_set<long int>& visited)
{
  for(auto &el : visited)
    {
      auto i = el/T;
      auto t = el%T;
      if (i == s)
        sbd.betweenness_exact[s][t] = sbd.betweenness_exact[s][t] - sbd.deltadot[s][t];
      sbd.betweenness_exact[i][t] = sbd.betweenness_exact[i][t] - sbd.deltasvvt[i][t];

    }
    //   for (int t = 0; t<T;t++) {
    //     sbd.betweenness_exact[s][t] = sbd.betweenness_exact[s][t] - sbd.deltadot[s][t];
    //   }

    // for (int i = 0; i < n; i++)
    //   {
    //     for (int t = 0; t<T;t++) {
    //       sbd.betweenness_exact[i][t] = sbd.betweenness_exact[i][t] - sbd.deltasvvt[i][t];
    //     }
    //   }
}

void print_pred_neighbour(Predecessor &G, const akt::Graph& g)
{
  int T = g.events.size();
  for (auto &it : G.ordered_neighb )
    {
      std::cout << "node  " << it.first/T << g.events[it.first%T] << "\n";
      for (auto &itt: it.second)
        {
          std::cout << "   time" << g.events[itt.first] << "\n";
          for(auto &ittt : itt.second)
            std::cout << "      node" << ittt << "\n";
        }
    }
}

std::vector<long int> optimalUpdateBetweenness(int s, const akt::Graph& g, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph&), bool (*cmp)(double, double), std::string walk_type, bool no_loop)
{
  std::pair<std::unordered_set<int>, std::unordered_set<int>> p;
  Predecessor G = Predecessor(g, sbd.pre, s);
  // if we are sure the cost function does not create loopos in the pred graph we can skip this phase
  if (!no_loop)
    {
      p = RemoveInfiniteFromPredecessor(s, G, sbd, cost, cmp, walk_type, g);
      int T = g.events.size();
      for(auto &elem : p.first)
        {
          sbd.sigmadot[elem/T][elem%T] = std::numeric_limits<double>::infinity(); 
        }
      p.first.insert(p.second.begin(), p.second.end());
      for(auto &e : p.first)
        {
          //      std::cout << "ici" << e/T << " " << e%T << "\n";
          sbd.sigma[e/T][e%T] = std::numeric_limits<double>::infinity();
          if (sbd.cur_best[e/T][e%T] == sbd.optimalNode[e/T])
            {
              sbd.totalSigma[e/T] += std::numeric_limits<double>::infinity();;
              sbd.totalSigmaT[e/T][e%T] = std::numeric_limits<double>::infinity();;
            }
        } 
    }

  sourcesSinksRemoveISolated(G,g);
  auto vis_sig = VolumePathAt(G, s, sbd, g);
  if(walk_type == "active")
      vis_sig = OptimalSigma(s, G, sbd, g, cost, p.first);
  else
    copySigmas(vis_sig,sbd,g);

  ComputeDeltaSvvt(G, s, sbd, g);
  if(walk_type == "active")
    CompleteDelta(G, sbd, g);
  //display_tot(sbd);
  PredecessorGraphToOrdered(G, g.events.size(), g.N(),walk_type);
  //  printf("fin ordered\n");
  //print_pred_neighbour(G, g);
  std::map<int,int> preced = BeforeNodes(G, g);
  //  display_tot(sbd);
  auto visited = GeneralContribution(g, G, s, sbd, preced, walk_type);
  UpdateBetweenness(sbd,  g.events.size(), visited);
  UpdateBetweenness_exact(sbd, g.events.size(),s, visited);
  return vis_sig;
}

// main algorithm for shortest paths variants in active case
std::vector<long int> optimalUpdateBetweennessBoost(int s, const akt::Graph& g, OptimalBetweennessData& sbd, std::string walk_type)
{
  Predecessor G = Predecessor(g, sbd.pre, s);
  //printPred(G, g);
  std::pair<std::unordered_set<int>, std::unordered_set<int>> p;
  sourcesSinksRemoveISolated(G,g);
  //auto vis_sig = VolumePathAt(G, s, sbd, g);
  auto vis_sig = OptimalSigmaBoost(s, G, sbd, g);
  ComputeDeltaSvvt(G, s, sbd, g);
  CompleteDelta(G, sbd, g);
  PredecessorGraphToOrdered(G, g.events.size(), g.N(),walk_type);
  std::map<int,int> preced;
  preced = BeforeNodes(G, g);
  std::unordered_set<long int> visited;
  //  display_tot(sbd);
  visited = GeneralContribution(g, G, s, sbd, preced, walk_type);
  UpdateBetweenness(sbd,  g.events.size(), visited);
  UpdateBetweenness_exact(sbd, g.events.size(),s, visited);
  return vis_sig;
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
    if (x < y)
      return true;
    return false;
  }
  bool larger(double x, double y)
  {
    if (x > y)
      return true;
    return false;
  }

  // Computes the betweenness measures
  std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double> > optimalBetweenness(const Graph& g, bool strict, std::string cost, std::string cmp, std::string walk_type, int numberNodes, bool bellman, bool bfs, bool no_loop)
  {
    double (*cost2)(Path*, int, const akt::Graph&);
    bool (*cmp2)(double, double);
    // printf("**********************************************************************\n");
    if (cost == "shortestfastest")
        cost2 = &co_sfp;
    else if(cost == "shortest")
      cost2 = &co_short;
    else if(cost == "foremost")
      cost2 = &co_first_arrival;
    else if(cost == "fastest")
      cost2 = &co_fastest;
    else if(cost == "shortestforemost")
      cost2 = &co_shortest_foremost;
    else if(cost == "shortestrestless")
      cost2 = &co_shortest_2restless;

    cmp2 = &larger;
    if (cmp == "le")
        cmp2 = &smaller;
    // std::cout << "lol" << g.N()<< g.T();
    auto sbd = OptimalBetweennessData(g);
    // Stores the betweenness values; initialize to 1 because of the formula for betweenness having a constant +1
    int s = 0;
    std::unordered_set<int> sampled;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, g.N()-1);
    int i = 0;

    if(bellman)
      {
        //run Bellman-Ford
        std::cout << "Bellman-Ford version! \n";
      }
    else if(bfs)
      {
        //run bfs
        std::cout << "BFS version! \n";
      }
    else
      {
        //run Dijkstra
        std::cout << "Dijkstra version! \n";
      }


    std::cout << "pred has no loops speed up ? " << no_loop << " *0 means no 1 means yes*\n";
    while(i < g.N() && i != numberNodes){
      if(numberNodes != -1)
        {
          s = distrib(gen);
          while(sampled.count(s) == 1)
            s = distrib(gen);
          sampled.insert(s);
        }
      //    for (int s = 0; s < g.N(); ++s) {
      printf("********************************* v %d / %d \n",s,g.N()-1);
      //display_tot(sbd);
      std::vector<long int> vis;
      if(bellman)
        {
          //run Bellman-Ford
          vis = BF(g, s, sbd, cost2, cmp2, walk_type);
        }
      else if(bfs)
        {
          //run bfs
          vis = BFS(g, s, sbd, cost2, cmp2, walk_type);
        }
      else
        {
          //run Dijkstra
          vis = dijkstra(g, s, sbd, cost2, cmp2, walk_type); 
        }

      // display_tot(sbd);
      auto vis2 = optimalUpdateBetweenness(s, g, sbd, cost2, cmp2, walk_type, no_loop);
      // display_tot(sbd);
      if (walk_type == "active")
        reinitializeHelperStructOptimal(g, s, sbd, vis2);
      else
        reinitializeHelperStructOptimal(g, s, sbd, vis);

      // display_tot(sbd);

      //      std::cout << "end parts "<< "\n" << std::flush;
      s++;
      i++;
    }
    totalBetweenness_compute(sbd, g.N(), g.events.size());
    //    std::cout << "end total_bet_comp "<< "\n" << std::flush;
    //display_tot(sbd);
    return {sbd.betweenness , sbd.betweenness_exact, sbd.totalBetweenness};

  }
  // computes the betweenness centrality on shortest paths variants
  std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>> , std::vector<double> > optimalBetweennessBoost(const Graph& g, bool strict, std::string cost, std::string walk_type, int numberNodes, int percentTime)
  {
    std::cout << "Boost version! \n";
    double untilTime;
    if(percentTime == -1)
      untilTime = std::numeric_limits<double>::infinity();
    else
      untilTime = g.minActualTime() + ((g.maxActualTime() - g.minActualTime())*percentTime)/100.0;
    std::cout << "looking until " << untilTime << " start time " << g.minActualTime() << " end time "<< g.maxActualTime()<<  "\n";
    int foremost = 0;
    double k;
    if(cost == "shortest")
      k = std::numeric_limits<double>::infinity();
    else if(cost == "shortestrestless")
      k = (g.maxActualTime() - g.minActualTime()) * 10 / 100;
    else if(cost == "shortestforemost")
      foremost = 1;
    std::cout << "REST = " << k << " "<<g.minActualTime() << " "<< g.maxActualTime()<<  "\n";
    auto sbd = OptimalBetweennessData(g);
    // Stores the betweenness values; initialize to 1 because of the formula for betweenness having a constant +1
    int s = 0;
    std::unordered_set<int> sampled;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, g.N()-1);
    int i = 0;
    while(i < g.N() && i != numberNodes){
      if(numberNodes != -1){
          s = distrib(gen);
          while(sampled.count(s) == 1)
            s = distrib(gen);
          sampled.insert(s);
        }
      optimalComputeDistancesSigmasBoost(g, strict, s, sbd,  walk_type, k, foremost, untilTime);
      std::vector<long int> vis2;
      if(walk_type == "active")
        {
          vis2 = optimalUpdateBetweennessBoost(s, g, sbd, walk_type);
          //display_tot(sbd);
          reinitializeHelperStructOptimal(g, s, sbd, vis2);
        }
      else
        {
          shortestEmptyStackUpdateBetweennessBoost(s, sbd, g, foremost);
          //display_tot(sbd);
          reinitializeHelperStructOptimalBoost(g, s, sbd);
          //display_tot(sbd);
        }
      s++;
      i++;
    }
    totalBetweenness_compute(sbd, g.N(), g.events.size());
    // display_tot(sbd);
    return {sbd.betweenness , sbd.betweenness_exact, sbd.totalBetweenness};
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
  std::vector<double> shortestBetweennessStatic(const akt::Graph& g)
  {
    auto res = std::vector<double>(g.N(), 0);
    for (int s = 0; s < g.N(); ++s) {
      auto stack = std::vector<int>();
      auto preds = std::vector<std::vector<int>>(g.N(), std::vector<int>());
      auto sigma = std::vector<int>(g.N(), 0);
      auto dists = std::vector<int>(g.N(), -1);
      std::unordered_set<int> nodes;
      sigma[s] = 1;
      dists[s] = 0;
      std::queue<int> q;
      q.push(s);
      while (!q.empty()) {
        auto v = q.front();
        q.pop();
        stack.push_back(v);
        nodes.clear();
        for (int t = g.minimalTimestep(); (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[v][t].nextTimestep) {
          for (auto w : g.adjacencyList()[v][t].neighbours) {
            if(nodes.count(w) == 0)
              nodes.insert(w);
              }
        }
        for(auto &w : nodes){
            //        for (auto w : (g.adjacencyList()[v][0]).neighbours) {
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
        //delta[w] = 1;
        stack.pop_back();
        for (auto v : preds[w])
          {
            auto tmp = static_cast<double>(sigma[v]) / static_cast<double>(sigma[w])* (1 + delta[w]);
            delta[v] += tmp;
            if (v != s)
              res[v] += tmp;
          }
      }
    }
    return res;
  }



} // end namespace akt

