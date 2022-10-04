#include <predecessor_graph.h>
#include <algorithm>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <numeric>
#include <map>
#include <tuple>

Predecessor* PredecessorGraph(const akt::Graph& g, std::vector<std::map<int, std::unordered_set<akt::VertexAppearance>>> pre, int node)
{
  int n = g.N();
  int T = g.maximalTimestep();
  // std::unordered_set<akt::VertexAppearance>> ens;
  // for (int k = 0; k < g.N(); i++)
  //   {
  //     map<int, std::unordered_set<akt::VertexAppearance>>::iterator key;
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
  Predecessor G;
  G.g = NetworKit::Graph(n*T, false, true, false);
  for (int k = 0; k < n; i++)
    {
      map<int, std::unordered_set<akt::VertexAppearance>>::iterator key;
      for(key=pre[k].begin(); key!=pre[k].end(); ++key){
        {
          for (const auto& elem: pre[i][key]) {
            v = elem.v;
            t = elem.time;
            if !(k == node  &&  elem.v != -1){
                if (v == node)
                  {
                    G.g.addEdge(node*g.maximalTimestep() + key,k*g.maximalTimestep() + key);
                  }
                else
                  {
                    G.g.addEdge(v*g.maximalTimestep() + t, k*g.maximalTimestep() + key);
                  }
              }
          }
        }
      }
    }
  return &G

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
      std::vector<node> stack = bfs.getNodesSortedByDistance();
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

std::pair<std::unordered_set<int>, std::unordered_set<int>> RemoveInfiniteFromPredecessor(int s, Predecessor& G, OptimalBetweennessData& sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), std::string walk_type)
{
  NetworKit::StronglyConnectedComponents scc(G.g);
  scc.run();
  NetworKit::Graph cond = condensationGraph(G.g, scc);
  std::unordered_set<int> inf_scc;
  const auto sizes = scc.getPartition().subsetSizes();
  const auto ids = scc.getPartition().getSubsetIds();
  const auto map_id = scc.getPartition().subsetSizeMap();
  map<int, int>::iterator it;
  for (it = map_id.begin(); it != map_id.end(); it++)
    {
      // we keep ids of scc with more than one temporal vertex
      if(it->second > 1)
        inf_scc.inserts(it->first);
    }
  std::unordered_set<int>::iterator itr;
  for (itr = inf_scc.begin(); itr != inf_scc.end(); itr++)
    {
      NetworKit::BFS bfs(cond, (*itr), false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      for(NetworKit::node i : stack)
          inf_scc.inserts(*i);
    }
  auto temp_inf = std::unordered_set<int>;
  for (it = inf_scc.begin(); it != inf_scc.end(); itt++)
    {
      const auto elem = std::set<NetowrKit::index> getMembers(*it);
      std::set<NetowrKit::index>::iterator itt;
      for (itt = elem.begin(); itr != elem.end(); itr++)
        {
          G.g.removeNode(*itt);
          temp_inf.inserts(*itt);
        }
    }
  for(int i = 0; i < n*T; i++)
    {
      if (G.g.hasNode(i))
        {
          if (G.g.degreeOut(i) == 0)
            {
              G.sinks.inserts(i);
              if( G.g.degreeIn(i) == 0)
                {
                  G.sources.inserts(i);
                  G.g.removeNode(i);
                }
            }
          else
            {
              if ( G.g.degreeIn(i) == 0)
                G.sources.inserts(i)
            }
        }
    }
  auto clos_inf = infinite_closure(s, G.g, temp_inf, cost, sbd, cmp, n);

  return {temp_inf, clos_inf};
}

unordered_set<akt::VertexAppearance> Infinite_closure(Predecessor& G, std::vector<int> &events, list<int> &events_rev, OptimalBetweennessData &sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), double n, std::unordered_set<int> &node_inf)
{
  int T = events[events.size() -1];
  unordered_set<akt::VertexAppearance> res;
  for (const auto& elem: node_inf) {
    int i = events_rev[elem.t] + 1;
    while(i < events.size() && G.g.hasNode(elem.v * T + elem.t) == false)
      {
        if (cmp(cost(sbd.opt_walk[elem.v][elem.t]), events[i], n) , sbd.cur_best[elem.v][events[i]])
          {
            res.inserts(akt::VertexAppearance{v, events[i]});
          }
        i++;
      }
  }
  return res;
}


void volumePathAtRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, int T)
{
  if(sbd.sigmadot[akt::VertexAppearance{e/T,e%T}] != 0)
    return;
  if(e/T == s)
    sbd.sigmadot[akt::VertexAppearance{e/T,e%T}] = 1;
  else
    {
      unsigned long long res = 0;
      G.g->forInEdgesOf(e, [&](node u, edgeweight w)
      {
        volumePathAtRec(s,u,G.g,sbd,T);
        if(u/T == s)
          res += 1;
        else
          res += sbd.sigmadot[u];
      });
      sigma[e] = res;
      if (sbd.cur_best[e/T][e%T] == sbd.optimalNode[e/T])
        {
          sbd.totalSigma[e/T] += res;
          sbd.totalSigmaT[e/T][e%T] = res;
        }
    }
  return;

}


void VolumePathAt(Predecessor& G, int s, set<NetworKit::index> sinks, OptimalBetweennessData &sbd)
{
  std::set<NetowrKit::index>::iterator itt;
  for (itt = G.sinks.begin(); itr != G.sinks.end(); itr++)
      volumePathAtRec(s,e,G.g,sbd);
}

void OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& s, double n, double (*cost)(Path, int, double))
{
  int T = s.maximalTimestep();
  for(int k = 0; k < s.N(); k++)
    {
      int pred = -1;
      for (int ev = 0; ev < s.events.size(); ev ++) {
          int t = s.events[ev];
          if (node == k)
            {
              if (sbd.sigmadot[akt::VertexAppearance{k, t}] > 0 )
                sbd.sigma[akt::VertexAppearance{k, t}] = sbd.sigmadot[akt::VertexAppearance{k, t}];
              else
                sbd.sigma[akt::VertexAppearance{k, t}] = 0;
            }
          else
            {
              if(pred == -1)
                {
                  if (G.g.hasNode(k*T + t))
                    {
                      sbd.sigma[akt::VertexAppearance{k, t}] = sbd.sigmadot[akt::VertexAppearance{k, t}];
                      pred = t;
                    }
                  else
                      sbd.sigma[akt::VertexAppearance{k, t}] = 0;
                }
              else
                {
                  if (G.g.hasNode(k*T + t))
                    {
                      if(sbd.cur_best[k][t] == cost(sbd.opt_walk[k][pred], t, n))
                        {
                          sbd.sigma[akt::VertexAppearance{k, t}] = sbd.sigma[akt::VertexAppearance{k, pred}] + sbd.sigmadot[akt::VertexAppearance{k, t}];
                          pred = t;
                        }
                      else
                        {
                          sbd.sigma[akt::VertexAppearance{k, t}] = sbd.sigmadot[akt::VertexAppearance{k, t}];
                          pred = t;
                        }
                    }
                  else
                    {
                      sbd.sigma[akt::VertexAppearance{k, t}] = sbd.sigma[akt::VertexAppearance{k, pred}];
                    }
                }
            }
          if (node_inf.count(akt::VertexAppearance{k,t}) == 1)
            {
              sbd.sigma[akt::VertexAppearance{k, t}] = std::numeric_limits<double>::infinity();
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

  G.g->forInEdgesOf(e, [&](node u, edgeweight w){computeDeltaRec(s,u,G.g,sbd,T);});
  deltasvvt[e/T][e%T] = sbd.totalSigmaT[e/T][e%T] / sbd.totalSigma[e/T];
  return;
}


void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd)
{
  std::set<NetowrKit::index>::iterator itt;
  for (itt = G.sinks.begin(); itr != G.sinks.end(); itr++)
    computeDeltaRec(s,e,G.g,sbd);
}


void PredecessorGraphToOrdered(Predecessor& G, int  T, std::map<int, int> ev_rev)
{
  //   std::map<int, std::map<int, std::vector<int> > > ordered_neighb;
  G.g.forEdges(
             [&](node x, node y, edgeid eid)
             {
               int v = x/T;
               int t = x%T;
               int w = y/T;
               int tp = y%T;
               if(G.g.ordered_neighb.count(x) == 1)
                 {
                   if (G.g.ordered_neighb[tp].count(tp) == 1)
                     G.g.ordered_neighb[x][ev_rev[tp]].push_back(w);
                   else
                     {
                       std::vector<int> m;
                       m.push_back(w);
                       G.g.ordered_neighb[x][ev_rev[tp]] = m;
                     }
                 }
               else
                 {
                   std::map<int, vector<int>> m;
                   std::vector<int> p;
                   p.push_back(w);
                   G.g.ordered_neighb[x][ev_rev[tp]] = m;
                 }


             });

}


void GeneralContribution(const akt::Graph& g, Predecessor G, int s, OptimalBetweennessData& sbd , std::map<VertexAppearance, VertexAppearance> &preced, std::string walk_type)
{
  //check if need to order
  for(auto star : G.sources) {
    int v = star/T;
    int t = star%T;
    unordered_set<VertexAppearance> visited;
    DeltaSvt(s, v, t, G.ordered_neighb, sbd, g, preced, walk_type, visited);
  }
}

void DeltaSvt(int s, int v, int t, std::map<int, std::map<int, std::vector<int> > >& l_nei, OptimalBetweennessData& sbd, const akt::Graph& g, std::map<VertexAppearance, VertexAppearance> &preced , std::string walk_type, unordered_set<VertexAppearance> visited)
{
  int T = g.maximalTimestep();
  if (visited.contains(VertexAppearance {v,t}))
    return;
  std::map<int, double> partial_sum;
  std::map<int, double> contrib_local;
  auto s = 0.0;
  std::map<int, std::map<int, std::vector<int> > >::reverse_iterator rit;
  for (rit=l_nei.rbegin(); rit!=l_nei.rend(); ++rit)
    {
      for(const auto& u: l_nei[v*T + t][rit->second])
        {
          int w = u/T;
          int tp = u%T;
          DeltaSvt(s, w, tp,  l_nei,  sbd, g, preced , walk_type, visited);
          s += (sbd.sigma[v][t]/sbd.sigma[w][tp])*sbd.deltadot[w][tp];
          partial_sum[rit->second] = s;
        }
    }
  sbd.deltadot[v][t] = s + sbd.deltasvvt[v][t];
}
