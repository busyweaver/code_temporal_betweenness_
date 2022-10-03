#include <algorithm>
#include<networkit/graph/Graph.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <numeric>
#include <map>
void predecessor_graph(const akt::Graph& g, std::vector<std::map<int, std::unordered_set<VertexAppearance>>> pre, int node)
{
  int n = g.N();
  int T = g.maximalTimestep();
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
  NetworKit::Graph G = NetworKit::Graph(n*T, false, true, false);
  for (int k = 0; k < n; i++)
    {
      map<int, std::unordered_set<VertexAppearance>>::iterator key;
      for(key=pre[k].begin(); key!=pre[k].end(); ++key){
        {
          for (const auto& elem: pre[i][key]) {
            v = elem.v;
            t = elem.time;
            if !(k == node  &&  elem.v != -1){
                if (v == node)
                  {
                    G.addEdge(node*g.maximalTimestep() + key,k*g.maximalTimestep() + key);
                  }
                else
                  {
                    G.addEdge(v*g.maximalTimestep() + t, k*g.maximalTimestep() + key);
                  }
              }
          }
        }
      }
    }

}

NetworKit::Graph condensationGraph(NetworKit::Graph& G, NetworKit::StronglyConnectedComponents scc)
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
      NetworKit::BFS bfs(G, *itr2, false, true);
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

void removeInfiniteFromPredecessor(int s, NetworKit::Graph& G, OptimalBetweennessData& sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), std::string walk_type, list[int] &events, list[int]& events_rev)
{
  NetworKit::StronglyConnectedComponents scc(G);
  scc.run();
  NetworKit::Graph cond = condensationGraph(G, scc);
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
        {
          inf_scc.inserts(*i);
        }
    }
  auto temp_inf = std::unordered_set<int>;
  for (it = inf_scc.begin(); it != inf_scc.end(); itt++)
    {
      const auto elem = std::set<NetowrKit::index> getMembers(*it);
      std::set<NetowrKit::index>::iterator itt;
      for (itt = elem.begin(); itr != elem.end(); itr++)
        {
          G.removeNode(*itt);
          temp_inf.inserts(*itt);
        }

    }
  for(int i = 0; i < n*T; i++)
    {
      if (G.hasNode(i))
        {
          if (G.degreeOut(i) == 0)
            {
              sinks.inserts(i);
              if( G.degreeIn(i) == 0)
                {
                  sources.inserts(i);
                  G.removeNode(i);
                }
            }
          else
            {
              if ( G.degreeIn(i) == 0)
                sources.inserts(i)
            }

        }
    }
  auto clos_inf = infinite_closure(s, G, temp_inf, cost, sbd, cmp, n);

  return {temp_inf, clos_inf, sources, sinks};
}

unordered_set<VertexAppearance> infinite_closure(NetworKit::Graph& G, list<int> &events, list<int> &events_rev, OptimalBetweennessData &sbd, double (*cost)(Path, int, double), bool (*cmp)(double, double), double n, std::unordered_set<int> &node_inf)
{
  int T = events[events.size() -1];
  unordered_set<VertexAppearance> res;
  for (const auto& elem: mySet) {
    int i = events_rev[elem.t] + 1;
    while(i < events.size() && G.hasNode(elem.v * T + elem.t) == false)
      {
        if (cmp(cost(sbd.opt_walk[elem.v][elem.t]), events[i], n) , sbd.cur_best[elem.v][events[i]])
          {
            res.inserts(VertexAppearance{v, events[i]});
          }
        i++;
      }
  }
  return res;
}


void volumePathAtRec(int s,int e,NetworKit::Graph& G,OptimalBetweennessData &sbd, int T)
{
  if(sbd.sigmadot[VertexAppearance{e/T,e%T}] == -1)
    return;
  if(e/T == s)
    sbd.sigmadot[VertexAppearance{e/T,e%T}] = 1;
  else
    {
      unsigned long long res = 0;
      G->forInEdgesOf(e, [&](node u, edgeweight w)
      {
        volumePathAtRec(s,u,G,sbd,T);
        if(u/T == s)
          res += 1;
        else
          res += sbd.sigmadot[u];
      });
      sigma[e] = res;
    }
  return;

}


void volumePathAt(NetworKit::Graph& G, int s, set<NetworKit::index> sinks, OptimalBetweennessData &sbd)
{
  std::set<NetowrKit::index>::iterator itt;
  for (itt = elem.begin(); itr != elem.end(); itr++)
    {
      volumePathAtRec(s,e,G,sbd)
    }
}

void optimal_sigma(int node, list[int] events, NetworKit::Graph &G, OptimalBetweennessData &sbd, const akt::Graph& s, double n, double (*cost)(Path, int, double))
{
  int T = G.maximalTimestep();
  for(int k = 0; k < s.N(); k++)
    {
      int pred = -1;
      for (int ev = 0; ev < events.size(); ev ++) {
          int t = events[ev];
          if (node == k)
            {
              if (sbd.sigmadot[VertexAppearance{k, t}] > 0 )
                sbd.sigma[VertexAppearance{k, t}] = sbd.sigmadot[VertexAppearance{k, t}];
              else
                sbd.sigma[VertexAppearance{k, t}] = 0;
            }
          else
            {
              if(pred == -1)
                {
                  if (G.hasNode(k*T + t))
                    {
                      sbd.sigma[VertexAppearance{k, t}] = sbd.sigmadot[VertexAppearance{k, t}];
                      pred = t;
                    }
                  else
                      sbd.sigma[VertexAppearance{k, t}] = 0;
                }
              else
                {
                  if (G.hasNode(k*T + t))
                    {
                      if(sbd.cur_best[k][t] == cost(sbd.opt_walk[k][pred], t, n))
                        {
                          sbd.sigma[VertexAppearance{k, t}] = sbd.sigma[VertexAppearance{k, pred}] + sbd.sigmadot[VertexAppearance{k, t}];
                          pred = t;
                        }
                      else
                        {
                          sbd.sigma[VertexAppearance{k, t}] = sbd.sigmadot[VertexAppearance{k, t}];
                          pred = t;
                        }
                    }
                  else
                    {
                      sbd.sigma[VertexAppearance{k, t}] = sbd.sigma[VertexAppearance{k, pred}];
                    }
                }
            }
          if (node_inf.count(VertexAppearance{k,t}) == 1)
            {
              sbd.sigma[VertexAppearance{k, t}] = std::numeric_limits<double>::infinity();
            }
        }

    }
}
