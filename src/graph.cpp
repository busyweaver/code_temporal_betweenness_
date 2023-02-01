#include "graph.h"

#include <sstream>
#include <unordered_map>
#include <tuple>
#include<iostream>
// Removes timesteps (i.e. "compresses" the time) where there are no temporal edges
void removeBoringTimesteps(akt::TemporalEdgeSet& edges, std::map<int, int> &events_rev)
{
    int timestepsRemoved = 0, lastTimestep = 0;
    // Handle first edge separately
    if (!edges.empty())
        timestepsRemoved = edges.cbegin()->when;
    
    // Flatten the set into an array to be able to easily modify the when-values
    auto flat = std::vector<akt::TemporalEdge>(edges.cbegin(), edges.cend());
    
    for (auto& te : flat) {

        // te.when -= timestepsRemoved;
        // if ((te.when - lastTimestep) > 1) {
        //     int newTime = lastTimestep + 1;
        //     timestepsRemoved += te.when - newTime;
        //     te.when = newTime;
        // }
        // lastTimestep = te.when;
      te.when = events_rev[te.when];
    }
    
    // Re-blow into a tree
    edges = akt::TemporalEdgeSet(flat.cbegin(), flat.cend(), &akt::temporalEdgeLessTimewise);
}

namespace akt {
    bool temporalEdgeLessTimewise(const akt::TemporalEdge& lhs, const akt::TemporalEdge& rhs)
    {
        return  std::tie(lhs.when, lhs.from, lhs.to) < std::tie(rhs.when, rhs.from, rhs.to);
        //return (lhs.when != rhs.when) ? (lhs.when < rhs.when)
        //    : ((lhs.from != rhs.from) ? (lhs.from < rhs.from)
        //    : (lhs.to < rhs.to));
    }

  std::pair<Graph, std::vector<std::string>> readReduceGraph(std::istream& is, bool directed, double eps)
    {
      std::set<int> events_set;
      std::vector<int> events;
      std::map<int, int> events_rev;

        TemporalEdgeSet edgesRead(&temporalEdgeLessTimewise);
        std::unordered_map<std::string, int> nodeIds;
        std::vector<std::string> reverseNodeIds;
        int noNodes = 0;
        for (std::string line; std::getline(is, line); ) {
            auto iss = std::istringstream{ line };
            //std::cout << "ligne "<< line;
            std::string from, to;
            int w;
            iss >> from >> to >> w;
            //std::cout << from << to << w;
            events_set.insert(w);
            if ((from.empty()) || (to.empty()))
                continue;
            if (nodeIds.count(from) < 1) {
                nodeIds[from] = noNodes++;
                reverseNodeIds.push_back(from);
            }
            if (nodeIds.count(to) < 1) {
                nodeIds[to] = noNodes++;
                reverseNodeIds.push_back(to);
            }
            int f = nodeIds[from];
            int t = nodeIds[to];

            if (f == t)
                continue;
            edgesRead.insert(TemporalEdge{ f, t, w });
            if (!directed)
                edgesRead.insert(TemporalEdge{ t, f, w });
        }
        int min = *(events_set.begin());
        int max = *(events_set.rbegin());
        if (eps > 0)
          {
            int j = min;
            while(j < max)
              {
                events_set.insert(j);
                j = j + eps;
              } 
          }
        events.assign(events_set.begin(), events_set.end());
        std::sort(events.begin(), events.end());

        for (int i = 0; i < events.size(); i++)
          events_rev[events[i]] = i;
        // printf("events now \n");

        removeBoringTimesteps(edgesRead, events_rev);

        return { Graph(noNodes, edgesRead.crbegin()->when, edgesRead, events, events_rev), reverseNodeIds };
    }
}
