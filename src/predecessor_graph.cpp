#include <algorithm>
#include "Snap.h"


void predecessor_graph(const akt::Graph& g, std::vector<std::vector<std::unordered_set<VertexAppearance>>> pre, int node)
{
  int n = g.N()
  PNGraph Graph = TNGraph::New();
  std::unordered_set<VertexAppearance>> ens;
  for (int i = 0; i < g.N(); i++)
    {
      for (int key = 0; key < g.maximalTimestep(); key ++)
        {
          for (const auto& elem: pre[i][key]) {
            v = elem.v;
            t = elem.time;
            if !(k == node  &&  elem.v != -1){
                if (v == node){
                  if (v == node)
                    {
                      if !(ens.contains(i*g.maximalTimestep() + key ) )
                        Graph->AddNode(VertexAppearance{k,key});

                    }
                }
              }
        }
    }
}
