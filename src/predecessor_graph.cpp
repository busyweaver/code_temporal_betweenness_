#include <algorithm>
#include "Snap.h"


void predecessor_graph(const akt::Graph& g, std::vector<std::map<int, std::unordered_set<VertexAppearance>>> pre, int node)
{
  int n = g.N()
  PNGraph Graph = TNGraph::New();
  std::unordered_set<VertexAppearance>> ens;
  for (int k = 0; k < g.N(); i++)
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
                      Graph->AddNode(k*g.maximalTimestep() + key);
                      Graph->AddNode(node*g.maximalTimestep() + key);
                      Graph->AddEdge(k*g.maximalTimestep() + key, node*g.maximalTimestep() + key);
                    }
                else
                  {
                    Graph->AddNode(k*g.maximalTimestep() + key);
                    Graph->AddNode(v*g.maximalTimestep() + t);
                    Graph->AddEdge(k*g.maximalTimestep() + key, v*g.maximalTimestep() + t);
                  }
                }
              }
        }
    }
}
}
