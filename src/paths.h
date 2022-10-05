#pragma once
#include "graph.h"
#include <limits>
#include <stdlib.h>
class Path {
private:
  std::vector<int>* nodes;
  std::vector<int>* times;
  int end;
public:
  Path(std::vector<int> * nod, std::vector<int> * time);
  bool is_empty();
  bool is_none();

  void add_link(int a,int b,int t);

  int arrival();
  int departure();
  int length();
  void display();

  Path clone();

};

double co_sfp(Path m, int t, const akt::Graph &g);

double co_short(Path m, int t, double cons);

double co_first_arrival(Path m, int t, double cons);
