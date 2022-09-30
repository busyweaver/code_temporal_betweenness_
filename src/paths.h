#include <limits>
#include <stdlib.h>
class Path {
private:
  std::vector<int>* nodes;
  std::vector<int>* times;
  int end;
public:
  Path(std::vector<int> * nod, std::vector<int> * time)
  {nodes = nod; times = time;
    end = (*times).size();}
  bool is_empty() {return (nodes != NULL) && ((*nodes).size() == 0);}
  bool is_none() {return nodes == NULL;}

  void add_link(int a,int b,int t)
  {
    if (nodes == NULL || times == NULL)
      exit(-1);
    if (nodes != NULL && (*nodes).size() == 0)
      {
        (*nodes).push_back(a);
        (*nodes).push_back(b);
        (*times).push_back(t);
      }
    else
      {
        if ((*nodes)[(*nodes).size() -1] == a)
          {
            (*nodes).push_back(b);
            (*times).push_back(t);
          }
      }
    end ++;

  }
  int arrival()
  {
    if(nodes == NULL || (*nodes).size() == 0)
      exit(-1);
    return (*times)[end -1];
  }
  int departure()
  {
    if(nodes == NULL || (*nodes).size() == 0)
      exit(-1);
    return (*times)[0];
  }
  int length()
  {
    if(nodes == NULL)
      exit(-1);
    return end;
  }
  void display()
  {
    for (int i = 0; i < (*times).size(); i++)
      {
        printf("%d - %d-> %d ", (*nodes)[i],(*times)[i],(*nodes)[i+1]);
      }

  }
  Path clone()
  {
    return Path(nodes, times);
  }


};

double co_sfp(Path m, int t, double cons)
{
  if (m.is_none())
    {
      return std::numeric_limits<double>::infinity();
    }
  if (m.is_empty())
    return 0.0;
  return cons * (t - m.departure()) + m.length();
}
double co_short(Path m, int t, double cons)
{
  if (m.is_none())
      return std::numeric_limits<double>::infinity();
  return m.length();
}

double co_first_arrival(Path m, int t, double cons)
{
  if (m.is_none())
      return std::numeric_limits<double>::infinity();
  if (m.is_empty())
    return 0.0;
  return m.arrival();
}

