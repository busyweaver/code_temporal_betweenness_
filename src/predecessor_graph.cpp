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

std::pair<std::map<long int,long int >,std::map<long int,long int >> numberNodes(const akt::Graph& gg, std::map<int, std::map<int,std::unordered_set<VertexAppearance>>> pre, int node)
{
  std::map<long int, long int> ma;
  std::map<long int, long int> starting;

  int actual = 0;
  int n = gg.N();
  int TT = gg.events.size();

  for(auto &el : pre)
    {
      auto k = el.first;
      for(auto &ell : pre[k])
        {
          auto key = ell.first;
          for (const auto& elem: pre[k][key]) {
            auto v = elem.v;
            auto t = elem.time;
            // std::cout << "node k " << k << " key " << key << " v " << v << " t " << t << "\n";
            if(!(k == node && t == -1))
              { 
                if (v == node)
                  {
                    // std::cout << node <<  gg.events[key] << "with "<< k <<  gg.events[key] << "\n";
                    if (starting.count(k) == 1)
                      {
                        if(starting[k] > key)
                          starting[k] = key;
                      }
                    else
                      starting[k] = key;
                    if (starting.count(node) == 1)
                      {
                        if(starting[node] > key)
                          starting[node] = key;
                      }
                    else
                      starting[node] = key;
                    if(!(ma.count(node*TT + key) == 1))
                      {
                        ma[node*TT + key] = actual;
                        actual ++;

                      }
                    if(!(ma.count(k*TT + key) == 1) )
                      {
                        ma[k*TT + key] = actual;
                        actual ++;
                      }

                  }
                else
                  {
                    //                    std::cout << "number " << k << " "<< key << " "<<starting.count(k)<< "\n"; 
                    if (starting.count(k) == 1)
                      {
                        if(starting[k] > key)
                          starting[k] = key;
                      }
                    else
                      starting[k] = key;
                    //                    std::cout << "number res =>" << k << " "<< key << " "<<starting[k]<< "\n"; 
                    // std::cout << v << gg.events[t] << "with "<< k << gg.events[key] << "\n";
                    if(!(ma.count(k*TT + key) == 1))
                      {
                        ma[k*TT + key] = actual;
                        actual ++;
                      }
                    if(!(ma.count(v*TT + t) == 1))
                      {
                        ma[v*TT + t] = actual;
                        actual ++;
                        // if (starting.count(v) == 1)
                        //   {
                        //     if(starting[v] > key)
                        //       starting[v] = key;
                        //   }
                      }

                  }
              }

          }
        }
    }
  // printf("nb events %d\n", TT);
  return {ma, starting};
}

Predecessor::Predecessor(const akt::Graph& gg, std::map<int, std::map<int,std::unordered_set<VertexAppearance>>> pre, int node)
{
  int actual = 0;
  int n = gg.N();
  int TT = gg.events.size();
  auto p = numberNodes(gg,pre,node);
  ma = p.first;
  starting_time = p.second;

  for (auto &e : ma)
    ma_inv[e.second] = e.first;

  g = NetworKit::Graph(ma.size(), false, true, false);
  for(auto &el : pre)
    {
      auto k = el.first;
        for(auto &ell : pre[k])
          {
            auto key = ell.first;
            for (const auto& elem: pre[k][key]) {
              auto v = elem.v;
              auto t = elem.time;
              // std::cout << "node k " << k << " key " << key << " v " << v << " t " << t << "\n";
              if(!(k == node && t == -1))
                {
                  if (v == node)
                    {
                      // std::cout << node <<  gg.events[key] << "with "<< k <<  gg.events[key] << "\n";

                      g.addEdge(ma[node*TT + key], ma[k*TT + key]);
                    }
                  else
                    {
                      // std::cout << v << gg.events[t] << "with "<< k << gg.events[key] << "\n";
                      g.addEdge(ma[v*TT + t], ma[k*TT + key]);
                    }
                }
            }
          }
    }
}

void printPred(Predecessor& G, const akt::Graph & g)
{
  auto T = g.events.size();
  G.g.forNodes(
               [&](NetworKit::node x_i)
               {
                 auto x = G.ma_inv[x_i];
                 auto v = x/T;
                 auto t = x%T;
                 std::cout << "node -> " << v << " " << t << "\n";
               });
  G.g.forNodes(
               [&](NetworKit::node v)
               {
                 std::cout << G.ma_inv[v]/T << " " <<  G.ma_inv[v]%T << "\n" << std::flush;
                 G.g.forNeighborsOf(v,
                                    [&](NetworKit::node x_i, NetworKit::edgeweight w)
                                    {
                                      std::cout << " -> " << G.ma_inv[x_i]/T << " " << G.ma_inv[x_i]%T << "\n" << std::flush;
                                    });
               });

}

NetworKit::Graph condensationGraph(Predecessor& G, NetworKit::StronglyConnectedComponents scc)
{
  //  std::cout  <<   "start conden \n" << std::flush;
  const auto partition = scc.getPartition();
  //  std::cout  <<   "1 \n" << std::flush;
  const auto par = scc.getPartition().getSubsets();
  //  std::cout  <<   "2 \n" << std::flush;
  const auto sizes = scc.getPartition().subsetSizes();

  const auto comp_sizes = scc.getComponentSizes();
  //  std::cout  <<   "3 \n" << std::flush;
  const auto ids = scc.getPartition().getSubsetIds();
  //std::set<index> getMembers(index s) const

  std::unordered_set<NetworKit::index> not_visited;
  for (auto &v : ids)
    not_visited.insert(v);

  // for (auto &v : comp_sizes)
  //   {
  //     if(v.second != 1)
  //       std::cout  << "sizes "<< v.second <<  " \n" << std::flush;
  //   }

  //  std::cout  <<   "init conden \n" << std::flush;
  NetworKit::Graph GG = NetworKit::Graph(sizes.size(), false, true, false);
  //  std::cout << "sizes visited" << not_visited.size() << " \n" << std::flush;
  while(not_visited.size() > 0)
    {
      //  std::cout  <<   "size not visired "<< not_visited.size() <<" \n" << std::flush;
      auto first_comp = not_visited.begin();
      auto set_comp =   partition.getMembers(*first_comp);
      //      std::cout  <<   "size actual comp "<< set_comp.size() <<" \n" << std::flush;
      auto first_elem = set_comp.begin();
      NetworKit::BFS bfs(G.g, *first_elem, false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      if(stack.size() > 1)
        //        std::cout << "   stack size " << stack.size() << " \n" << std::flush;
      for(NetworKit::node i : stack)
        {
          //  std::cout  <<   "   stack " << i << "\n" << std::flush;
          if( partition.subsetOf(i) !=  partition.subsetOf(*first_elem))
            {
              not_visited.erase(partition.subsetOf(i));
              GG.addEdge(partition.subsetOf(*first_elem),partition.subsetOf(i));
            }
        }
      not_visited.erase(partition.subsetOf(*first_elem));
    }

  // std::set<std::set<NetworKit::index>>::iterator itr;
  
  // for (itr = par.begin(); itr != par.end(); itr++)
  //   {
  //     std::set<NetworKit::index>::iterator itr2;
  //     itr2 = (*itr).begin();
  //     //      std::cout  <<   "start bfs on graph " << *itr2 << "\n" << std::flush;
  //     NetworKit::BFS bfs(G.g, *itr2, false, true);
  //     bfs.run();
  //     //      std::cout  <<   "end bfs on graph \n" << std::flush;
  //     std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
  //     for(NetworKit::node i : stack)
  //       {
  //         //  std::cout  <<   "   stack " << i << "\n" << std::flush;
  //         if( partition.subsetOf(i) !=  partition.subsetOf(*itr2))
  //           {
  //             GG.addEdge(partition.subsetOf(*itr2),partition.subsetOf(i));
  //           }
  //       }
      //std::cout << "i = " << i << std::endl;

  //  }
  return GG;
}


std::unordered_set<int> Infinite_closure(Predecessor& G, OptimalBetweennessData &sbd, double (*cost)(Path*, int, const akt::Graph &), bool (*cmp)(double, double), std::unordered_set<int> &node_inf, const akt::Graph &g)
{
  // std::cout << "init infinite closure : ev size" << g.events.size() << "\n";
  // for (auto &elem: g.events) {
  //   std::cout << "event " << elem << "\n" << std::flush;
  // }
  auto T = g.events.size();
  // std::cout << "T" << T << " \n"<< std::flush;
  std::unordered_set<int> res;
  for (auto &elem: node_inf) {
    // std::cout << elem <<"node inf : "<< elem/T << " "<< elem%T  <<   "\n" << std::flush;

    auto i = g.events_rev.at(g.events[elem%T]) + 1;
    // std::cout << " i " << i <<  " elem/T + g.events[i] "<< elem/T + g.events[i] <<"\n" << std::flush;
    while(i < g.events.size() && G.g.hasNode(elem/T + g.events[i]) == false)
      {
        // std::cout << " elem/T, elem%T  "<< elem/T  << " " <<  elem%T <<"\n" << std::flush;
        // should not it be sbd.cur_best[elem/T][elem%T] ?!
        if (cmp(cost(sbd.opt_walk[elem/T][elem%T], g.events[i], g) , sbd.cur_best[elem/T][g.events[i]]))
          {
            res.insert(elem/T + g.events[i]);
          }
        i++;
      }
  }
  // std::cout << "end infinite closure \n" << std::flush;

  return res;
}


void sourcesSinksRemoveISolated(Predecessor& G, const akt::Graph & g)
{
  G.g.forNodes(
               [&](NetworKit::node i)
               {
                 if (G.g.degreeOut(i) == 0)
                   {
                     if( G.g.degreeIn(i) == 0)
                       {
                         std::cout  <<   "removing \n" << std::flush;  
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
               });
  //  std::cout  <<   "end sources and sinks \n" << std::flush;  
}

std::pair<std::unordered_set<int>, std::unordered_set<int>> RemoveInfiniteFromPredecessor(int s, Predecessor& G, OptimalBetweennessData& sbd, double (*cost)(Path*, int, const akt::Graph &), bool (*cmp)(double, double), std::string walk_type, const akt::Graph & g)
{
  NetworKit::StronglyConnectedComponents scc(G.g);
  scc.run();
  //  std::cout  <<   "end connected comp \n" << std::flush;
  NetworKit::Graph cond = condensationGraph(G, scc);
  //  std::cout  <<   "end conden\n" << std::flush;
  std::unordered_set<int> inf_scc;
  const auto  partition = scc.getPartition();
  //  std::cout  <<   "end partition\n" << std::flush;
  const auto sizes = scc.getPartition().subsetSizes();
  //  std::cout  <<   "end sizes\n" << std::flush;
  const auto ids = scc.getPartition().getSubsetIds();
  //  std::cout  <<   "end subset ids \n" << std::flush;
  const auto map_id = scc.getPartition().subsetSizeMap();
  //  std::cout  <<   "end subset sizemap \n" << std::flush;
  std::map<NetworKit::index, NetworKit::count>::const_iterator it;
  // std::cout << "init remove \n";
  for (it = map_id.begin(); it != map_id.end(); it++)
    {
      // we keep ids of scc with more than one temporal vertex
      if(it->second > 1)
        inf_scc.insert(it->first);
    }
  //  std::cout << "scc of size > 1   " << inf_scc.size() <<  "\n";
  std::unordered_set<int>::iterator itr;
  for (itr = inf_scc.begin(); itr != inf_scc.end(); itr++)
    {
      NetworKit::BFS bfs(cond, (*itr), false, true);
      bfs.run();
      std::vector<NetworKit::node> stack = bfs.getNodesSortedByDistance();
      for(auto &i : stack)
          inf_scc.insert(i);
    }
  //std::cout  <<   "end bfs on comp \n" << std::flush;
  // std::cout << "trans scc \n";
  std::unordered_set<int > temp_inf;
  for (auto &itr : inf_scc)
    {
      // std::cout << "comp scc /**/ \n"; 
      const auto elem = partition.getMembers(itr);
      long int T = g.events.size();
      for (auto &itt : elem)
        {
          //          std::cout << "inf node"<<  G.ma_inv[itt]/T << " " << G.ma_inv[itt]%T<< "\n"; 
          G.g.removeNode(itt);
          temp_inf.insert(G.ma_inv[itt]);
        }
    }
  //std::cout  <<   "end remove infinite without closure \n" << std::flush;

  // std::cout << "remove isolated?\n";
  std::unordered_set<int> clos_inf;
  if(walk_type == "active")
    clos_inf = Infinite_closure(G, sbd, cost, cmp, temp_inf, g);

  //  for (auto &itt :temp_inf)
      // std::cout << "SHOULD inf node"<<  itt << "\n";
  //for (auto &itt :clos_inf)
    // std::cout << "SHOULD2 inf node"<<  itt << "\n";

  // std::cout << "left before close "<<  G.g.numberOfNodes() << "\n";
  //  std::cout << "size close inf "<<  clos_inf.size() << "\n";
  //  std::cout << "size temp inf "<<  temp_inf.size() << "\n";
  for (auto &itt :clos_inf)
    {
      if (G.g.hasNode(G.ma[itt]))
        G.g.removeNode(G.ma[itt]);
    }
  //  std::cout  <<   "end remove all infinite \n" << std::flush;
  //printPred(G, g);
  sourcesSinksRemoveISolated(G, g);
  //printPred(G, g);
  return {temp_inf, clos_inf};
}

void copySigmas(std::unordered_set<long int> &visited,OptimalBetweennessData &sbd, const akt::Graph &g)
{
  int T = g.events.size();
  for(auto &e : visited)
    sbd.sigma[e/T][e%T] = sbd.sigmadot[e/T][e%T];

}


void volumePathAtRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, const akt::Graph &g, std::unordered_set<long int> &visited)
{

  int T = g.events.size();
  //  std::cout << "start rec sigma original node " << s << "node " << e/T << " "<< e%T << " \n";
  if(sbd.sigmadot[e/T][e%T] != 0)
    return;
  if(G.g.degreeIn(G.ma[e]) == 0)
    {
      //maybe only source node should be put to 1 by adding this source as an argument to this function
      // std::cout << "ici\n";
      sbd.sigmadot[e/T][e%T] = 1;
      //      sbd.sigma[e/T][e%T] = 1;
      visited.insert(e);
    }

  else
    {
      // std::cout << "from  "<< e/T << " "<<e%T << "\n";
      double res = 0;
      G.g.forInEdgesOf(G.ma[e], [&](NetworKit::node u_i, NetworKit::edgeweight w)
      {
        auto u = G.ma_inv[u_i];

        volumePathAtRec(s,u,G,sbd,g, visited);

        // if(u/T == s)
        //   res += 1;
        // else
          res += sbd.sigmadot[u/T][u%T];

      });
      // printf("fin\n");
      sbd.sigmadot[e/T][e%T] = res;
      //      sbd.sigma[e/T][e%T] = res;
      visited.insert(e);
      // std::cout << "volume dot " <<e/T << " " << e%T << " = " << res << "\n";

      if (sbd.cur_best[e/T][e%T] == sbd.optimalNode[e/T])
        {
          sbd.totalSigma[e/T] += res;
          sbd.totalSigmaT[e/T][e%T] = res;
        }
    }
}


std::unordered_set<long int> VolumePathAt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph &g)
{
  std::unordered_set<long int> visited;
  for (auto &itt : G.sinks)
    {
      //      printf("ssssssssssssssiiiiiiiiiiiiiinnnnnnnnnnnks %ld %ld\n",G.ma_inv[itt]/g.events.size(), G.ma_inv[itt]%g.events.size());
      volumePathAtRec(s,G.ma_inv[itt],G,sbd,g, visited);
    }
  return visited;
}

std::unordered_set<long int> OptimalSigma(int node, Predecessor &G, OptimalBetweennessData &sbd, const akt::Graph& g,  double (*cost)(Path*, int, const akt::Graph& g), std::unordered_set<int> node_inf)
{ //start function
  //  printf("nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn");
  std::unordered_set<long int> visited;
  int T = g.events.size();
  for (auto &elem : G.ma)
    {
      //      std::cout << "G -> "<< "T " << T << "elem "<< elem.first << " "<< elem.first/T << " " << elem.first%T << "\n";
    }
  //  std::cout << "1 1 -> " << " " <<  G.ma.count(T +1) <<"\n";
  for(int k = 0; k < g.N(); k++)
    {

      int pred = G.starting_time[k];
      for (int t = G.starting_time[k]; t < T ;t ++)
        {
          visited.insert(k*T + t);
          //          std::cout << k << " " << t << "\n";
          // int pred = -1;
          // for (int t = 0; t < T ;t ++) {
          if (node == k)
            {
              if (sbd.sigmadot[k][t] > 0 )
                sbd.sigma[k][t] = sbd.sigmadot[k][t];
              else
                sbd.sigma[k][t] = 0;
            }
          else
            {
              if (G.ma.count(k*T +t) == 1 &&  G.g.hasNode(G.ma[k*T + t]))
                {
                  //std::cout << "lol \n";
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
          if (node_inf.count(k*T + t) == 1)
            {
              sbd.sigma[k][t] = std::numeric_limits<double>::infinity();
            }

        }
    }
  return visited;
}


void CompleteDelta(Predecessor& G, OptimalBetweennessData &sbd, const akt::Graph& g)
{
  int T = g.events.size();
  for(int k = 0; k < g.N(); k++)
    {
      //for (int t = 0; t < T ;t ++) {
      for (long int t = G.starting_time[k]; t < T ;t ++) {
        if (!G.g.hasNode(k*T + t))
          {
            if (sbd.totalSigma[k] != 0)
              {
                sbd.deltasvvt[k][t] = sbd.totalSigmaT[k][t] / sbd.totalSigma[k];
                sbd.deltadot[k][t] = sbd.totalSigmaT[k][t] / sbd.totalSigma[k];
              }

          }
      }
    }
}

void computeDeltaRec(int s,int e,Predecessor& G,OptimalBetweennessData &sbd, int T, std::unordered_set<int>&visited)
{
  if(visited.count(e) == 1)
    return;
  if(e/T == s)
    return;// sbd.deltasvvt[e/T][e%T] = 0;
  else if (sbd.totalSigma[e/T] == std::numeric_limits<double>::infinity())
    return; //sbd.deltasvvt[e/T][e%T] = 0;
    else
      {
        sbd.deltasvvt[e/T][e%T] = sbd.totalSigmaT[e/T][e%T] / sbd.totalSigma[e/T];
        sbd.deltadot[e/T][e%T] = sbd.totalSigmaT[e/T][e%T] / sbd.totalSigma[e/T];
      }

  visited.insert(e);

  G.g.forInEdgesOf(G.ma[e], [&](NetworKit::node u, NetworKit::edgeweight w){computeDeltaRec(s,G.ma_inv[u],G,sbd,T, visited);});

  return;
}


void ComputeDeltaSvvt(Predecessor& G, int s, OptimalBetweennessData &sbd, const akt::Graph& g)
{
  std::unordered_set<int> visited;
  std::set<NetworKit::index>::iterator itt;
  for (auto &itt : G.sinks)
    computeDeltaRec(s,G.ma_inv[itt],G,sbd, g.events.size(), visited);
}


void PredecessorGraphToOrdered(Predecessor& G, int  T, int n, std::string walk_type)
{
  //   std::map<int, std::map<int, std::vector<int> > > ordered_neighb;
  G.g.forEdges(
               [&](NetworKit::node x_i, NetworKit::node y_i, NetworKit::edgeid eid)
             {
               auto x = G.ma_inv[x_i];
               auto y = G.ma_inv[y_i];
               int v = x/T;
               int t = x%T;
               int w = y/T;
               int tp = y%T;
               G.ordered_neighb[x][tp].push_back(w);

             });

  for (auto &el : G.ordered_neighb)
    {
      auto v = el.first;
      for (auto &tmp : G.ordered_neighb[v])
        {
          G.times_ord[v].push_back(tmp.first);
        }
      if(walk_type == "active")
        std::sort(G.times_ord[v].begin(), G.times_ord[v].end(), std::greater<long int>()); 
    }
}


std::map<int,int> BeforeNodes(Predecessor& G, const akt::Graph& g)
{
  int T = g.events.size();
  std::map<int,int> res;
  std::map<int,std::vector<int>> d;
  G.g.forNodes(
               [&](NetworKit::node x_i)
               {
                 auto x = G.ma_inv[x_i];
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

void IntermediaryNodes(int vt, int vtp, std::map<int,int>& before,OptimalBetweennessData& sbd, const akt::Graph& g, double s, int pred_time, std::unordered_set<long int>& visited)
{
  auto T = g.events.size();
  // std::cout  << "vt " << vt/T << " " << vt%T << "vtp " << vtp/T << " "<<vtp%T << "\n";
  auto v = vtp / T;
  auto tpp = vtp % T;
  auto t = vt % T;
  // std::cout << "before val : " << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "before[v*T + tpp ]" << before[v*T + tpp ]/T << " "<<before[v*T + tpp ]%T <<"\n" << std::flush;

  while (tpp >= pred_time && before[v*T + tpp ]%T == t && visited.count(v*T + tpp) == 0)
    {
      // std::cout << "*********** change val inter val" << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "\n" << std::flush;
      sbd.deltadot[v][tpp] = s + sbd.deltasvvt[v][tpp];
      tpp = tpp - 1;
      visited.insert(v*T + tpp);
      // std::cout << "before val" << "v "<< v << " tpp "<< tpp << " t " << t<< "pred_time " << pred_time << "\n" << std::flush;
    }
}

void DeltaSvt(int v, Predecessor& G, OptimalBetweennessData& sbd, const akt::Graph& g, std::map<int, int> &preced , std::string walk_type, std::unordered_set<long int>& visited)
{
  // for(auto &elem : visited)
  //   printf("**************************************************************///////\n");
  std::map<long int, std::map<long int, std::vector<long int> > >& l_nei = G.ordered_neighb;
  int T = g.events.size();
  if (visited.count(v) == 1)
    {
      // std::cout << "stop immediately" <<   "\n" << std::flush;;
      return;
    }
  //  std::cout << "new"  << v/T << " " << v%T << "\n";
  //  std::map<int, double> partial_sum;
  std::map<int, double> contrib_local;
  auto s = 0.0;

  for(int j = 0;j<G.times_ord[v].size();j++)
    //  for (auto &tp : times_ord)
    {
      int pred_time;
      auto tp = G.times_ord[v][j];
      if(j<G.times_ord[v].size()-1)
        pred_time = G.times_ord[v][j+1]+1;
      else
        pred_time = v%T+1;

      // std::cout << "check increasing time " << tp;
      for (auto &w : l_nei[v][tp])
        {
          //          std::cout << "next" << v/T << " " << v%T << " -> "  << w << " " << tp << "\n" << std::flush;
          // char tmp;
          // scanf("%c",&tmp);

          DeltaSvt(w*T + tp,  G,  sbd, g, preced , walk_type, visited);
          if(sbd.sigma[w][tp] == 0)
            {
              printf("gros probleme %d %d -> %ld %ld dot %lg sig %lg",v/T,v%T,w,tp,sbd.sigmadot[w][tp], sbd.sigma[w][tp]);
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
              // printf("ev : %d", ev);
              if(ev >= 0)
                IntermediaryNodes(v,(v/T)*T + ev,preced,sbd,g,s,pred_time, visited);
            }
    }
  sbd.deltadot[v/T][v%T] += s;// + sbd.deltasvvt[v/T][v%T];
  visited.insert(v);
  // std::cout << "end "  << v/T << " " << v%T << " s "<< s << "\n";
}


std::unordered_set<long int> GeneralContribution(const akt::Graph& g, Predecessor& G, int s, OptimalBetweennessData& sbd , std::map<int, int> &preced, std::string walk_type)
{
  int T = g.events.size();
  std::unordered_set<long int> visited;
  //check if need to order
  for(auto star : G.sources) {
    //    DeltaSvt(star, G.ordered_neighb, sbd, g, preced, walk_type, visited);
    DeltaSvt(G.ma_inv[star], G, sbd, g, preced, walk_type, visited);
  }

  return visited;
}

