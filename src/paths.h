#pragma once
#include "graph.h"
#include <limits>
#include <stdlib.h>
class Path {
public:
  int node;
  int time;
  int depar;
  int length;
  Path* parent;
// private:
//   std::vector<int>* nodes;
//   std::vector<int>* times;
//   int end;
// public:
//   Path(std::vector<int> * nod, std::vector<int> * time);

  Path(Path* p, int node, int time, int dep, int length);
  Path(Path* p, int node, int time);

  bool is_empty();
  bool is_none();
  int arrival();
  int departure();
  int getLength();


  Path clone();

};
void display(Path*);

double co_sfp(Path* m, int t, const akt::Graph &g);

double co_short(Path m, int t, double cons);

double co_first_arrival(Path m, int t, double cons);
