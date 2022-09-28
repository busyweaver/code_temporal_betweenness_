#pragma once

#include <utility>
#include <vector>

#include "graph.h"

namespace akt {
	// Returns the shortest betweenness and shortest foremost betweenness (with strictness depending on the parameter)
	std::pair<std::vector<double>, std::vector<double>> shortestBetweenness(const akt::Graph& g, bool strict);

	// Returns the strict prefix-foremost betwenness
	std::vector<double>	prefixForemostBetweenness(const akt::Graph& g);

	// Returns the shortest betweenness for a static graph (i.e. all edges have timestamp 0) using Brandes' algorithm
	std::vector<double> shortestBetweennessStatic(const akt::Graph& g);

} // end namespace akt