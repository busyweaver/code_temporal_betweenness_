#include "predecessor_graph.h"
#include <algorithm>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/graph/Graph.hpp>
#include <numeric>
#include <map>
#include <tuple>
//#include<vector>
Predecessor::Predecessor()
{}

Predecessor::Predecessor(const akt::Graph& gg, std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre, int node)
{
  
  int n = gg.N();
  int T = gg.events[gg.events.size()-1]+1;
  int TT = gg.events.size();
  // std::unordered_set<VertexAppearance>> ens;
  // for (int k = 0; k < g.N(); i++)
  //   {
  //     map<int, std::unordered_set<VertexAppearance>>::iterator key;
  //     for(key=pre[k].begin(); key!=pre[k].end(); ++key){
  //       {
  //         for (const auto& elem: pre[i][key]) {
  //           v = elem.v;
  //           t = elem.time;
  //           if !(k == node  &&  elem.v != -1){
  //               if (v == node)
  //                 {
  //                   ens.inserts(k*g.maximalTimestep() + key);
  //                   ens.inserts(node*g.maximalTimestep() + key);
  //                 }
  //               else
  //                 {
  //                   ens.inserts(k*g.maximalTimestep() + key);
  //                   ens.inserts(v*g.maximalTimestep() + t);
  //                 }
  //             }
  //         }
  //       }
  //     }
  //   }
  for (auto &ev : gg.events)
    std::cout << ev;
  printf("size graph TT %d n %d T %d\n", TT,n,T);
  g = NetworKit::Graph(n*(T), false, true, false);
  for (int k = 0; k < n; k++)
    {
      for(int key = 0; key < TT; key++){
        {
          for (const auto& elem: pre[k][key]) {
            auto v = elem.v;
            auto t = elem.time;
            std::cout << "elem " << v << " t " << t << " k " << k << " key " << key <<"events[key]" << gg.events[key] <<" starting " << node <<"\n";
            if(!(k == node && t == -1))
              { 
                if (v == node)
                  {
                    std::cout << node*T + (gg.events[key]) << node*T + (gg.events[key]) << "\n";
                    g.addEdge(node*T + (gg.events[key]), k*T + (gg.events[key]));
                  }
                else
                  {
                    std::cout << v*T + (t) <<  k*T + (gg.events[key]) << "\n";
                    g.addEdge(v*T + (t), k*T + (gg.events[key]));
                  }
              }
          }
        }
      }
    }
}

NetworKit::Graph condensationGraph(Predecessor& G, NetworKit::StronglyConnectedComponents scc)
{

  const auto partition = scc.getPartition();
  const auto par = scc.getPartition().getSubsets();
  const auto sizes = scc.getPartition().subsetSizes();
  const auto ids = scc.getPartition().getSubsetIds();
  NetworKit::Graph GG = NetworKit::Graph(sizes.size(), false, true, false);

  std::set<std::set<NetworKit::index>>::iterator itr;
  for (itr = par.begin(); itr != par.end(); itr++)
    {
      std::set<NetworKit::index>::iterator itr2;
      itr2 = (*itr).begin();
      NetworKit::BFS bfs(G.g, *itr2, false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      for(NetworKit::node i : stack)
        {
          if( partition.subsetOf(i) !=  partition.subsetOf(*itr2))
            {
              GG.addEdge(partition.subsetOf(*itr2),partition.subsetOf(i));
            }
        }
      //std::cout << "i = " << i << std::endl;

    }
  return GG;
}

std::unordered_set<int> Infinite_closure(Predecessor& G, OptimalBetweennessData &sbd, double (*cost)(Path*, int, const akt::Graph &), bool (*cmp)(double, double), std::unordered_set<int> &node_inf, const akt::Graph &g)
{
  int T = g.events[g.events.size() -1];
  std::unordered_set<int> res;
  for (auto &elem: node_inf) {
    int i = g.events_rev.at(elem%T) + 1;
    while(i < g.events.size() && G.g.hasNode(elem/T + g.events[i]) == false)
      {
        if (cmp(cost(sbd.opt_walk[elem/T][elem%T], g.events[i], g) , sbd.cur_best[elem/T][g.events[i]]))
          {
            res.insert(elem/T + g.events[i]);
          }
        i++;
      }
  }
  return res;
}



std::pair<std::unordered_set<int>, std::unordered_set<int>> RemoveInfiniteFromPredecessor(int s, Predecessor& G, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph &), bool (*cmp)(double, double), std::string walk_type, const akt::Graph & g)
{
  NetworKit::StronglyConnectedComponents scc(G.g);
  scc.run();
  NetworKit::Graph cond = condensationGraph(G, scc);
  std::unordered_set<int> inf_scc;
  const auto  partition = scc.getPartition();
  const auto sizes = scc.getPartition().subsetSizes();
  const auto ids = scc.getPartition().getSubsetIds();
  const auto map_id = scc.getPartition().subsetSizeMap();
  std::map<NetworKit::index, NetworKit::count>::const_iterator it;
  for (it = map_id.begin(); it != map_id.end(); it++)
    {
      // we keep ids of scc with more than one temporal vertex
      if(it->second > 1)
        inf_scc.insert(it->first);
    }
  std::unordered_set<int>::iterator itr;
  for (itr = inf_scc.begin(); itr != inf_scc.end(); itr++)
    {
      NetworKit::BFS bfs(cond, (*itr), false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      for(auto &i : stack)
          inf_scc.insert(i);
    }
  std::unordered_set<int > temp_inf;
  for (auto &itr : inf_scc)
    {
      const auto elem = partition.getMembers(itr);
      for (auto &itt : elem)
        {
          G.g.removeNode(itt);
          temp_inf.insert(itt);
        }
    }
  int T = g.events[g.events.size()-1]+1;
  for(int i = 0; i < g.N()*T; i++)
    {
      if (G.g.hasNode(i))
        {
          if (G.g.degreeOut(i) == 0)
            {
              if( G.g.degreeIn(i) == 0)
                {
                  G.g.removeNode(i);
                }
              else
                G.sinks.insert(i);
            }
          else
            {
              if ( G.g.degreeIn(i) == 0)
                G.sources.insert(i);
            }
        }
    }
  auto clos_inf = Infinite_closure(G, sbd, cost, cmp, temp_inf, g);

  return {temp_inf, clos_inf};
}


void volumePathAtRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, const akt::Graph &g)
{
  int T = g.events[g.events.size() -1]+1;
  if(sbd.sigmadot[e/T][g.events_rev.at(e%T)] != 0)
    return;
  if(e/T == s)
    sbd.sigmadot[e/T][g.events_rev.at(e%T)] = 1;
  else
    {
      unsigned long long res = 0;
      G.g.forInEdgesOf(e, [&](NetworKit::node u, NetworKit::edgeweight w)
      {
        volumePathAtRec(s,u,G,sbd,g);
        if(u/T == s)
          res += 1;
        else
          res += sbd.sigmadot[u/T][g.events_rev.at(u%T)];
      });
      sbd.sigmadot[e/T][g.events_rev.at(e%T)] = res;
      if (sbd.cur_best[e/T][g.events_rev.at(e%T)] == sbd.optimalNode[e/T])
        {
          sbd.totalSigma[e/T] += res;
          sbd.totalSigmaT[e/T][g.events_rev.at(e%T)] = res;
        }
    }
  return;

}


void VolumePathAt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph &g)
{
  printf("alo\n");
  for (auto &itt : G.sinks)
    {
      printf("sinks %d",itt);
      volumePathAtRec(s,itt,G,sbd,g);
    }

}

void OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& g,  double (*cost)(Path*, int, const akt::Graph& g), std::unordered_set<int> node_inf)
{ //start function
  int T = g.events[g.events.size()-1]+1;
  for(int k = 0; k < g.N(); k++)
    {
      int pred = -1;
      for (int ev = 0; ev < g.events.size(); ev ++) {
          int t = g.events[ev];
          if (node == k)
            {
              if (sbd.sigmadot[k][t] > 0 )
                sbd.sigma[k][t] = sbd.sigmadot[k][t];
              else
                sbd.sigma[k][t] = 0;
            }
          else
            {
              if(pred == -1)
                {
                  if (G.g.hasNode(k*T + t))
                    {
                      sbd.sigma[k][t] = sbd.sigmadot[k][t];
                      pred = t;
                    }
                  else
                      sbd.sigma[k][t] = 0;
                }
              else
                {
                  if (G.g.hasNode(k*T + t))
                    {
                      if(sbd.cur_best[k][t] == cost(sbd.opt_walk[k][pred], t, g))
                        {
                          sbd.sigma[k][t] = sbd.sigma[k][pred] + sbd.sigmadot[k][t];
                          pred = t;
                        }
                      else
                        {
                          sbd.sigma[k][t] = sbd.sigmadot[k][t];
                          pred = t;
                        }
                    }
                  else
                    {
                      sbd.sigma[k][t] = sbd.sigma[k][pred];
                    }
                }
            }
          if (node_inf.count(k*T + t) == 1)
            {
              sbd.sigma[k][t] = std::numeric_limits<double>::infinity();
            }
        }

    }
}

void computeDeltaRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, int T)
{
  if(sbd.deltasvvt[e/T][e%T] != 0)
    return;
  if(e/T == s)
    sbd.deltasvvt[e/T][e%T] = 0;
  if (sbd.totalSigma[e/T] = std::numeric_limits<double>::infinity())
    sbd.deltasvvt[e/T][e%T] = 0;

  G.g.forInEdgesOf(e, [&](NetworKit::node u, NetworKit::edgeweight w){computeDeltaRec(s,u,G,sbd,T);});
  sbd.deltasvvt[e/T][e%T] = sbd.totalSigmaT[e/T][e%T] / sbd.totalSigma[e/T];
  return;
}


void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph& g)
{
  std::set<NetworKit::index>::iterator itt;
  for (auto &itt : G.sinks)
    computeDeltaRec(s,itt,G,sbd, g.maximalTimestep());
}


void PredecessorGraphToOrdered(Predecessor& G, int  T, std::map<int, int> ev_rev)
{
  //   std::map<int, std::map<int, std::vector<int> > > ordered_neighb;
  G.g.forEdges(
               [&](NetworKit::node x, NetworKit::node y, NetworKit::edgeid eid)
             {
               int v = x/T;
               int t = x%T;
               int w = y/T;
               int tp = y%T;
               if(G.ordered_neighb.count(x) == 1)
                 {
                   if (G.ordered_neighb[tp].count(tp) == 1)
                     G.ordered_neighb[x][ev_rev[tp]].push_back(w);
                   else
                     {
                       std::vector<int> m;
                       m.push_back(w);
                       G.ordered_neighb[x][ev_rev[tp]] = m;
                     }
                 }
               else
                 {
                   // std::map<int, std::vector<int>> m;
                   // std::vector<int> p;
                   // p.push_back(w);
                   G.ordered_neighb[x][ev_rev[tp]].push_back(w);
                 }


             });

}
void DeltaSvt(int v, int t, std::map<int, std::map<int, std::vector<int> > >& l_nei, OptimalBetweennessData& sbd, const akt::Graph& g, std::map<VertexAppearance, VertexAppearance> &preced , std::string walk_type, std::unordered_set<VertexAppearance> visited)
{
  int T = g.maximalTimestep();
  if (visited.count(VertexAppearance {v,t}) == 1)
    return;
  std::map<int, double> partial_sum;
  std::map<int, double> contrib_local;
  auto s = 0.0;
  std::map<int, std::vector<int> >::reverse_iterator rit;
  for (rit=l_nei[v*T + t].rbegin(); rit!=l_nei[v*T + t].rend(); ++rit)
    {
      for(const auto &u: l_nei[v*T + t][(rit->first)])
        {
          int w = u;
          int tp = (rit->first);
          DeltaSvt(w, tp,  l_nei,  sbd, g, preced , walk_type, visited);
          s += (sbd.sigma[v][t]/sbd.sigma[w][tp])*sbd.deltadot[w][tp];
          partial_sum[rit->first] = s;
        }
    }
  sbd.deltadot[v][t] = s + sbd.deltasvvt[v][t];
}


void GeneralContribution(const akt::Graph& g, Predecessor G, int s, OptimalBetweennessData& sbd , std::map<VertexAppearance, VertexAppearance> &preced, std::string walk_type)
{
  int T = g.maximalTimestep();
  //check if need to order
  for(auto star : G.sources) {
    int v = star/T;
    int t = star%T;
    std::unordered_set<VertexAppearance> visited;
    DeltaSvt(v, t, G.ordered_neighb, sbd, g, preced, walk_type, visited);
  }
}

