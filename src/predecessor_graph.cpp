#include "predecessor_graph.h"
#include <algorithm>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/graph/Graph.hpp>
#include <numeric>
#include <map>
#include <tuple>
#include <ostream>      // std::flush

//#include<vector>
Predecessor::Predecessor()
{}

Predecessor::Predecessor(const akt::Graph& gg, std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre, int node)
{
  int n = gg.N();
  int TT = gg.events.size();
  for (auto &ev : gg.events)
    std::cout << ev;
  printf("size graph TT %d n %d\n", TT,n);
  g = NetworKit::Graph(n*(TT), false, true, false);
  for (int k = 0; k < n; k++)
    {
      for(int key = 0; key < TT; key++){
        {
          for (const auto& elem: pre[k][key]) {
            auto v = elem.v;
            auto t = elem.time;
            std::cout << "node k " << k << " key " << key << " v " << v << " t " << t << "\n";
            if(!(k == node && t == -1))
              { 
                if (v == node)
                  {
                    std::cout << node <<  gg.events[key] << "with "<< k <<  gg.events[key] << "\n";
                    g.addEdge(node*TT + key, k*TT + key);
                  }
                else
                  {
                    std::cout << v << gg.events[t] << "with "<< k << gg.events[key] << "\n";
                    g.addEdge(v*TT + (t), k*TT + key);
                  }
              }
          }
        }
      }
    }
  std::cout << "fin pred \n";
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
  std::cout << "init infinite closure : ev size" << g.events.size() << "\n";
  // for (auto &elem: g.events) {
  //   std::cout << "event " << elem << "\n" << std::flush;
  // }
  auto T = g.events.size();
  std::cout << "T" << T << " \n"<< std::flush;
  std::unordered_set<int> res;
  for (auto &elem: node_inf) {
    std::cout << elem <<"node inf : "<< elem/T << " "<< elem%T  <<   "\n" << std::flush;

    auto i = g.events_rev.at(g.events[elem%T]) + 1;
    std::cout << " i " << i <<  " elem/T + g.events[i] "<< elem/T + g.events[i] <<"\n" << std::flush;
    while(i < g.events.size() && G.g.hasNode(elem/T + g.events[i]) == false)
      {
        std::cout << " elem/T, elem%T  "<< elem/T  << " " <<  elem%T <<"\n" << std::flush;
        // should not it be sbd.cur_best[elem/T][elem%T] ?!
        if (cmp(cost(sbd.opt_walk[elem/T][elem%T], g.events[i], g) , sbd.cur_best[elem/T][g.events[i]]))
          {
            res.insert(elem/T + g.events[i]);
          }
        i++;
      }
  }
  std::cout << "end infinite closure \n" << std::flush;
  for (auto &er: res) {
    std::cout << "res " <<  er <<  "\n" << std::flush; 
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
  std::cout << "init remove \n";
  for (it = map_id.begin(); it != map_id.end(); it++)
    {
      // we keep ids of scc with more than one temporal vertex
      if(it->second > 1)
        inf_scc.insert(it->first);
    }
  std::cout << "scc of size > 1 \n";
  std::unordered_set<int>::iterator itr;
  for (itr = inf_scc.begin(); itr != inf_scc.end(); itr++)
    {
      NetworKit::BFS bfs(cond, (*itr), false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      for(auto &i : stack)
          inf_scc.insert(i);
    }
  std::cout << "trans scc \n";
  std::unordered_set<int > temp_inf;
  for (auto &itr : inf_scc)
    {
      std::cout << "comp scc /**/ \n"; 
      const auto elem = partition.getMembers(itr);
      for (auto &itt : elem)
        {
          std::cout << "inf node"<<  itt << "\n"; 
          G.g.removeNode(itt);
          temp_inf.insert(itt);
        }
    }
  std::cout << "remove nodes\n";
  int T = g.events.size();
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
  std::cout << "remove isolated?\n";
  auto clos_inf = Infinite_closure(G, sbd, cost, cmp, temp_inf, g);
  for (auto &itt :temp_inf)
      std::cout << "SHOULD inf node"<<  itt << "\n";
  for (auto &itt :clos_inf)
    std::cout << "SHOULD2 inf node"<<  itt << "\n";

  std::cout << "left before close "<<  G.g.numberOfNodes() << "\n";
  for (auto &itt :clos_inf)
    if (G.g.hasNode(itt))
       G.g.removeNode(itt);


  return {temp_inf, clos_inf};
}


void volumePathAtRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, const akt::Graph &g)
{

  int T = g.events.size();
  std::cout << "start rec sigma original node " << s << "node " << e/T << " "<< e%T << " \n";
  if(sbd.sigmadot[e/T][e%T] != 0)
    return;
  if(G.g.degreeIn(e) == 0)
    {
      //maybe only source node should be put to 1 by adding this source as an argument to this function
      std::cout << "ici\n";
      sbd.sigmadot[e/T][e%T] = 1;
    }

  else
    {
      std::cout << "from  "<< e/T << " "<<e%T << "\n";
      unsigned long long res = 0;
      G.g.forInEdgesOf(e, [&](NetworKit::node u, NetworKit::edgeweight w)
      {
        std::cout << "to "<< u/T << " "<< u%T << "\n";
        volumePathAtRec(s,u,G,sbd,g);
        if(u/T == s)
          res += 1;
        else
          res += sbd.sigmadot[u/T][u%T];
      });
      printf("fin\n");
      sbd.sigmadot[e/T][e%T] = res;
      std::cout << "volume dot " <<e/T << " " << e%T << " = " << res << "\n";

      if (sbd.cur_best[e/T][e%T] == sbd.optimalNode[e/T])
        {
          sbd.totalSigma[e/T] += res;
          sbd.totalSigmaT[e/T][e%T] = res;
        }
    }
  return;

}


void VolumePathAt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph &g)
{
  // printf("alo\n");
  // G.g.forNodes(
  //              [&](NetworKit::node x)
  //              {
  //                std::cout << "pred graph nodes : " << x/g.events.size() << " " << x%g.events.size() << " indegree "<<  G.g.degreeIn(x) << " outdegree " << G.g.degreeOut(x) << "\n" ;
  //              });
  // for (auto &itt : G.sinks)
  //     printf("sinks %ld %ld",itt/g.events.size(), itt%g.events.size());
  for (auto &itt : G.sinks)
    {
      printf("ssssssssssssssiiiiiiiiiiiiiinnnnnnnnnnnks %ld %ld",itt/g.events.size(), itt%g.events.size());
      volumePathAtRec(s,itt,G,sbd,g);
    }

}

void OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& g,  double (*cost)(Path*, int, const akt::Graph& g), std::unordered_set<int> node_inf)
{ //start function
  int T = g.events.size();
  for(int k = 0; k < g.N(); k++)
    {
      int pred = -1;
      for (int t = 0; t < T ;t ++) {
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
                      if(sbd.cur_best[k][t] == cost(sbd.opt_walk[k][pred], g.events[t], g))
                        {

                          sbd.sigma[k][t] = sbd.sigma[k][pred] + sbd.sigmadot[k][t];
                          pred = t;
                          if (sbd.cur_best[k][t] == sbd.optimalNode[k])
                            sbd.totalSigmaT[k][t] = sbd.sigma[k][t];
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
                      if (cost(sbd.opt_walk[k][pred], g.events[t], g) == sbd.optimalNode[k])
                        sbd.totalSigmaT[k][t] = sbd.sigma[k][t];
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


void CompleteDelta(Predecessor& G, OptimalBetweennessData &sbd, const akt::Graph& g)
{
  int T = g.events.size();
  for(int k = 0; k < g.N(); k++)
    {
      for (int t = 0; t < T ;t ++) {
        if (!G.g.hasNode(k*T + t))
          {
            if (sbd.totalSigma[k] != 0)
              sbd.deltasvvt[k][t] = sbd.totalSigmaT[k][t] / sbd.totalSigma[k];
          }
      }
    }
}

void computeDeltaRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, int T, std::unordered_set<int>&visited)
{
  if(visited.count(e) == 1)
    return;
  if(e/T == s)
    sbd.deltasvvt[e/T][e%T] = 0;
  else if (sbd.totalSigma[e/T] == std::numeric_limits<double>::infinity())
    sbd.deltasvvt[e/T][e%T] = 0;
  else
    sbd.deltasvvt[e/T][e%T] = sbd.totalSigmaT[e/T][e%T] / sbd.totalSigma[e/T];
  visited.insert(e);

  G.g.forInEdgesOf(e, [&](NetworKit::node u, NetworKit::edgeweight w){computeDeltaRec(s,u,G,sbd,T, visited);});

  return;
}


void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph& g)
{
  std::unordered_set<int> visited;
  std::set<NetworKit::index>::iterator itt;
  for (auto &itt : G.sinks)
    computeDeltaRec(s,itt,G,sbd, g.events.size(), visited);
}


void PredecessorGraphToOrdered(Predecessor& G, int  T)
{
  //   std::map<int, std::map<int, std::vector<int> > > ordered_neighb;
  G.g.forEdges(
               [&](NetworKit::node x, NetworKit::node y, NetworKit::edgeid eid)
             {
               int v = x/T;
               int t = x%T;
               int w = y/T;
               int tp = y%T;
               G.ordered_neighb[x][tp].push_back(w);

             });

}


std::map<int,int> BeforeNodes(Predecessor& G, const akt::Graph& g)
{
  int T = g.events.size();
  std::map<int,int> res;
  std::map<int,std::vector<int>> d;
  G.g.forNodes(
               [&](NetworKit::node x)
               {
                 auto v = x/T;
                 auto t = x%T;
                 d[v].push_back(t);
                 res[x] = x;
               });

  for (auto& [key, dv] : d)
    {
      std::sort(dv.begin(), dv.end());
      for (int i =0; i < dv.size() - 1; i++)
        {
          auto j = dv[i] + 1;
          while (j < dv[i+1])
            {
              res[key*T + j] = res[key*T + dv[i]];
              j = j + 1;
            }
        }
      auto j = dv[dv.size() - 1];
      while(j < g.events.size())
        {
          res[key*T + j] = res[key*T + dv[dv.size() - 1]];
          j = j+1;
        }
    }
  return res;
}

void IntermediaryNodes(int vt, int vtp, std::map<int,int> before,OptimalBetweennessData& sbd, const akt::Graph& g, double s, int pred_time)
{
  auto T = g.events.size();
  std::cout  << "vt " << vt/T << " " << vt%T << "vtp " << vtp/T << " "<<vtp%T << "\n";
  auto v = vtp / T;
  auto tpp = vtp % T;
  auto t = vt % T;
  std::cout << "before val : " << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "before[v*T + tpp ]" << before[v*T + tpp ]/T << " "<<before[v*T + tpp ]%T <<"\n" << std::flush;

  while (tpp >= pred_time && before[v*T + tpp ]%T == t)
    {
      std::cout << "*********** change val inter val" << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "\n" << std::flush;
      sbd.deltadot[v][tpp] = s + sbd.deltasvvt[v][tpp];
      tpp = tpp - 1;
      std::cout << "before val" << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "\n" << std::flush;
    }
}

void DeltaSvt(int v, Predecessor& G, OptimalBetweennessData& sbd, const akt::Graph& g, std::map<int, int> &preced , std::string walk_type, std::unordered_set<int>& visited)
{
  // for(auto &elem : visited)
  //   printf("**************************************************************///////\n");
  std::map<int, std::map<int, std::vector<int> > > l_nei = G.ordered_neighb;
  int T = g.events.size();
  std::cout << "new"  << v/T << " " << v%T << "\n";
  if (visited.count(v) == 1)
    {
      std::cout << "stop immediately" <<   "\n" << std::flush;;
      return;
    }
  //  std::map<int, double> partial_sum;
  std::map<int, double> contrib_local;
  auto s = 0.0;
  std::vector<int> times_ord;
  for (auto &tmp : l_nei[v])
    {
      times_ord.push_back(tmp.first);
    }
  std::sort(times_ord.begin(), times_ord.end(), std::greater<int>());
  for(int j = 0;j<times_ord.size();j++)
    //  for (auto &tp : times_ord)
    {
      int pred_time;
      auto tp = times_ord[j];
      if(j<times_ord.size()-1)
        pred_time = times_ord[j+1]+1;
      else
        pred_time = v%T+1;

      std::cout << "check increasing time " << tp;
      for (auto &w : l_nei[v][tp])
        {
          std::cout << "next" << v/T << " " << v%T << " -> "  << w << " " << tp << "\n" << std::flush;
          // char tmp;
          // scanf("%c",&tmp);

          DeltaSvt(w*T + tp,  G,  sbd, g, preced , walk_type, visited);
          if(sbd.sigma[w][tp] == 0)
            {
              printf("gros probleme %d %d -> %d %d",v/T,v%T,w,tp);
              exit(-1);
            }
          s += (sbd.sigma[v/T][v%T]/sbd.sigma[w][tp])*sbd.deltadot[w][tp];
          //          partial_sum[tp] = s;
        }
      if(walk_type == "active")
        {
          int ev = tp;
          if(G.g.hasNode((v/T)*T + tp))
            ev = tp -1;
          printf("ev : %d", ev);
          if(ev >= 0)
            IntermediaryNodes(v,(v/T)*T + ev,preced,sbd,g,s,pred_time); 
        }
    }
  sbd.deltadot[v/T][v%T] = s + sbd.deltasvvt[v/T][v%T];
  visited.insert(v);
  std::cout << "end "  << v/T << " " << v%T << " s "<< s << "\n";
}


void GeneralContribution(const akt::Graph& g, Predecessor& G, int s, OptimalBetweennessData& sbd , std::map<int, int> &preced, std::string walk_type)
{
  int T = g.events.size();
  //check if need to order
  for(auto star : G.sources) {
    std::unordered_set<int> visited;
    //    DeltaSvt(star, G.ordered_neighb, sbd, g, preced, walk_type, visited);
    DeltaSvt(star, G, sbd, g, preced, walk_type, visited);
  }
}

