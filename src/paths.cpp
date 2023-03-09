#include "paths.h"
#include<iostream>
Path::Path(Path* p, int nod, int tim, int dep, int len)
{
  this->parent = p;
  this->node = nod;
  this->time = tim;
  this->depar = dep;
  this->length = len;
}

Path::Path(Path* p,int b,int t)
{

  if (p == nullptr)
    {
      //printf("parent null\n");
      this->node = b;
      this->time = t;
      this->parent = p;
      this->depar = t;
      this->length = 0;
      //printf("size %d \n", this->length);
      //fflush(stdout);
    }
  else
    {
      //printf("parent ok\n");
      this->node = b;
      this->time = t;
      this->parent = p;
      // if(p -> parent == nullptr)
      //   this->depar = t;
      this -> depar = parent -> depar;
      this-> length = this->parent->length +1;

    }
  //printf("ici4\n");
  //fflush(stdout);
}
// Path::Path(std::vector<int> * nod, std::vector<int> * time)
// {this->nodes = nod; this->times = time;
//   if (this->nodes == nullptr)
//     end = -1;
//   else
//     end = (*this->times).size();
// }

bool Path::is_empty() {return this->node==-1;}
bool Path::is_none() {return this == nullptr;}



int Path::arrival()
{
  return this->time;
}

int Path::departure()
{
  return this->depar;
}

int Path::getLength()
{
  return this->length;
}

void display(Path*m)
{
  Path* tmp = m;
  if(tmp==nullptr)
    {
      std::cout << "nein empty \n";
      return;
    }
  int i =0;
  printf("\nPath\n");
  while(tmp->parent != nullptr && i<10)
    {
      printf("(%d <-(%d)- %d) ", tmp->node,tmp->time,tmp->parent->node);
      tmp = tmp->parent;
      i++;
    }
  std::cout << "\nInfo path depar" << (*m).departure() << "arrival " << (*m).arrival() << "length " << (*m).getLength() << " \n";
  printf("\nEnd Path\n");
}


double co_sfp(Path* m, int t, const akt::Graph &g)
{
  if (m == nullptr)
    {
      return std::numeric_limits<double>::infinity();
    }
  if ((*m).is_empty())
    return 0.0;
  //std::cout << "t" << t << "(*m).departure())" << (*m).departure() << "(*m).getLength()" << (*m).getLength() << "\n";
  return g.N() * (t - (*m).departure()) + (*m).getLength();
}


double co_short(Path *m, int t,  const akt::Graph &g)
{
  if (m == nullptr)
    return std::numeric_limits<double>::infinity();
  return (*m).getLength();
}

//for active types it is not prefix-optimal
double co_first_arrival(Path *m, int t,  const akt::Graph &g)
{
  if (m == nullptr)
    return std::numeric_limits<double>::infinity();
  if ((*m).is_empty())
    return 0.0;
  return (*m).arrival();
}

//for active types it is not prefix-optimal
double co_shortest_foremost(Path *m, int t,  const akt::Graph &g)
{
  if (m == nullptr)
    return std::numeric_limits<double>::infinity();
  if ((*m).is_empty())
    return 0.0;
  return g.N() * (*m).arrival() +   (*m).getLength();
}

double co_shortest_2restless(Path *m, int t,  const akt::Graph &g)
{
  double res;
  //  display(m);
  if (m == nullptr)
    res = std::numeric_limits<double>::infinity();
  else if ((*m).getLength() == 1 || (*m).getLength() == 0)
    res = 1;
  else
    {
      if(m->time - m->parent->time >2)
        res = std::numeric_limits<double>::infinity();
      else
        res = (*m).getLength();
    }
  //  std::cout << res << " res \n";
  return res;
}
