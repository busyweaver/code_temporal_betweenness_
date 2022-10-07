#include "paths.h"

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
  printf("path %p b %d t %d\n",(void*)p,b,t);
  printf("debut\n");
  fflush(stdout);
  // if (this->nodes == nullptr || this->times == nullptr)
  //   {
  //     printf("probleme\n");
  //     fflush(stdout);
  //     exit(-1);
  //   }

  if (p == nullptr)
    {
      printf("parent null\n");
      this->node = b;
      this->time = t;
      this->parent = p;
      this->depar = t;
      this->length = 0;
      printf("size %d \n", this->length);
      fflush(stdout); 
    }
  else
    {
      printf("parent ok\n");
      this->node = b;
      this->time = t;
      this->parent = p;
 

      this-> length = this->parent->length +1;
      this->depar = p->depar;
    }
  printf("ici4\n");
  fflush(stdout);
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
  int i =0;
  printf("%p,%p\n",(void*)m, (void*)m->parent);
  while(tmp->parent != nullptr && i<10)
    {
      printf("(%d <-(%d)- %d) ", tmp->node,tmp->time,tmp->parent->node);
      tmp = tmp->parent;
      i ++;
    }
}


double co_sfp(Path* m, int t, const akt::Graph &g)
{
  if (m == nullptr)
    {
      return std::numeric_limits<double>::infinity();
    }
  if ((*m).is_empty())
    return 0.0;
  return g.N() * (t - (*m).departure()) + (*m).getLength();
}


double co_short(Path m, int t, double cons)
{
  if (m.is_none())
    return std::numeric_limits<double>::infinity();
  return m.getLength();
}


double co_first_arrival(Path m, int t, double cons)
{
  if (m.is_none())
    return std::numeric_limits<double>::infinity();
  if (m.is_empty())
    return 0.0;
  return m.arrival();
}
