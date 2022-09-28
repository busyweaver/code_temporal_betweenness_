#include <limits>
#include <stdlib.h>
template struct Path {
private:
  void * nodes;
  void * times;
public:
  Path() {nodes = NULL;
    times = NULL;}
  Path(node, time)
  {nodes = node; times = time;}
  bool is_empty() {return nodes.size() == 0;}
  bool is_none() {return nodes == NULL;}

  void add_link(a int ,b int ,t int)
  {
    if (nodes != NULL && nodes.size() == 0)
      {
        nodes.push_back(a);
        nodes.push_back(b);
        times.push_back(t);
      }
    else
      {
        if (nodes[nodes.size() -1] == a)
          {
            nodes.push_back(b);
            times.push_back(t);
          }
      }

  }
  double co_sfp(int t, double cons)
  {
    if (nodes == NULL)
      {
        return std::numeric_limits<double>::infinity();
      }
    if (nodes.size() == 0)
        return 0.0
    return cons * (t - times[0]) + times.size();
  }
  double co_short(int t, double cons)
  {
    if (nodes == NULL)
      {
        return std::numeric_limits<double>::infinity();
      }
    return times.size();
   }

  double co_first_arrival(int t, double cons)
  {
    if (nodes == NULL)
      {
        return std::numeric_limits<double>::infinity();
      }
    return times[times.size() -1];
  }

  void display()
  {
    for (int i = 0; i < times.size(); i++)
      {
        printf("%d - %d-> %d ", nodes[i],times[i],nodes[i+1]);
      }

  }


};
