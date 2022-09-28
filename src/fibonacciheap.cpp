/* Implementation from https://rosettacode.org/wiki/Fibonacci_heap#C++ */
#include "fibonacciheap.h"

FibonacciHeap<V>::FibonacciHeap() {
  heap=_empty();
}
FibonacciHeap<V>::virtual ~FibonacciHeap() {
  if(heap) {
    _deleteAll(heap);
  }
}
FibonacciHeap<V>::node<V>* insert(V value) {
  node<V>* ret=_singleton(value);
  heap=_merge(heap,ret);
  return ret;
}
FibonacciHeap<V>::void merge(FibonacciHeap& other) {
  heap=_merge(heap,other.heap);
  other.heap=_empty();
}

FibonacciHeap<V>::bool isEmpty() {
  return heap==NULL;
}

FibonacciHeap<V>::V getMinimum() {
  return heap->value;
}

FibonacciHeap<V>::V removeMinimum() {
  node<V>* old=heap;
  heap=_removeMinimum(heap);
  V ret=old->value;
  delete old;
  return ret;
}

FibonacciHeap<V>::void decreaseKey(node<V>* n,V value) {
  heap=_decreaseKey(heap,n,value);
}

FibonacciHeap<V>::node<V>* find(V value) {
  return _find(heap,value);
}

FibonacciHeap<V>::void display() const {	// function code adapted from GO code just below C++
  node<V>* p = heap;
  if (p == NULL) {
    cout << "The Heap is Empty" << endl;
    return;
  }
  cout << "The root nodes of Heap are: " << endl;
  _display_tree(heap, "");
  cout << endl;
}



/*Private part*/

FibonacciHeap<V>::node<V>* _empty() {
  return NULL;
}

FibonacciHeap<V>::node<V>* _singleton(V value) {
  node<V>* n=new node<V>;
  n->value=value;
  n->prev=n->next=n;
  n->degree=0;
  n->marked=false;
  n->child=NULL;
  n->parent=NULL;
  return n;
}

FibonacciHeap<V>::node<V>* _merge(node<V>* a,node<V>* b) {
  if(a==NULL)return b;
  if(b==NULL)return a;
  if(a->value>b->value) {
    node<V>* temp=a;
    a=b;
    b=temp;
  }
  node<V>* an=a->next;
  node<V>* bp=b->prev;
  a->next=b;
  b->prev=a;
  an->prev=bp;
  bp->next=an;
  return a;
}

FibonacciHeap<V>::void _deleteAll(node<V>* n) {
  if(n!=NULL) {
    node<V>* c=n;
    do {
      node<V>* d=c;
      c=c->next;
      _deleteAll(d->child);
      delete d;
    } while(c!=n);
  }
}

FibonacciHeap<V>::void _addChild(node<V>* parent,node<V>* child) {
  child->prev=child->next=child;
  child->parent=parent;
  parent->degree++;
  parent->child=_merge(parent->child,child);
}

FibonacciHeap<V>::void _unMarkAndUnParentAll(node<V>* n) {
  if(n==NULL)return;
  node<V>* c=n;
  do {
    c->marked=false;
    c->parent=NULL;
    c=c->next;
  }while(c!=n);
}

FibonacciHeap<V>::node<V>* _removeMinimum(node<V>* n) {
  _unMarkAndUnParentAll(n->child);
  if(n->next==n) {
    n=n->child;
  } else {
    n->next->prev=n->prev;
    n->prev->next=n->next;
    n=_merge(n->next,n->child);
  }
  if(n==NULL)return n;
  node<V>* trees[64]={NULL};

  while(true) {
    if(trees[n->degree]!=NULL) {
      node<V>* t=trees[n->degree];
      if(t==n)break;
      trees[n->degree]=NULL;
      t->prev->next=t->next;
      t->next->prev=t->prev;
      if(n->value<t->value) {
        _addChild(n,t);
      } else {
        if(n->next==n) {
          t->next=t->prev=t;
        } else {
          n->prev->next=t;
          n->next->prev=t;
          t->next=n->next;
          t->prev=n->prev;
        }
        _addChild(t,n);
        n=t;
      }
      continue;
    } else {
      trees[n->degree]=n;
    }
    n=n->next;
  }
  node<V>* min=n;
  do {
    if(n->value<min->value)min=n;
    n=n->next;
  } while(n!=n);
  return min;
}

FibonacciHeap<V>::node<V>* _cut(node<V>* heap,node<V>* n) {
  if(n->next==n) {
    n->parent->child=NULL;
  } else {
    n->next->prev=n->prev;
    n->prev->next=n->next;
    n->parent->child=n->next;
  }
  n->next=n->prev=n;
  n->marked=false;
  n->parent->degree--;
  return _merge(heap,n);
}

FibonacciHeap<V>::node<V>* _decreaseKey(node<V>* heap,node<V>* n,V value) {
  if(n->value<value)return heap;
  n->value=value;
  node<V>* parent = n->parent;
  if(parent != nullptr && n->value < parent->value) {
    heap=_cut(heap,n);
    n->parent=NULL;
    while(parent!=NULL && parent->marked) {
      heap=_cut(heap,parent);
      n=parent;
      parent=n->parent;
      n->parent=NULL;
    }
    if(parent!=NULL && parent->parent!=NULL)parent->marked=true;
    if (n->value < heap->value)heap = n;
  }
  return heap;
}

FibonacciHeap<V>::node<V>* _find(node<V>* heap,V value) {
  node<V>* n=heap;
  if(n==NULL)return NULL;
  do {
    if(n->value==value)return n;
    node<V>* ret=_find(n->child,value);
    if(ret)return ret;
    n=n->next;
  }while(n!=heap);
  return NULL;
}

FibonacciHeap<V>::void _display_tree(node<V>* n, string pre) const {
  string pc = "│  ";
  node<V>* x = n;
  do {
    if (x->next != n) {
      cout << pre << "├─";
    } else {
      cout << pre << "└─";
      pc = "   ";
    }
    if (x->child == nullptr) {
      cout << "─" << x->value << endl;
    } else {
      cout << "┐" << x->value << endl;
      _display_tree(x->child, pre + pc);
    }
    //		cout << endl;
    x = x->next;
  } while (x != n);
}

