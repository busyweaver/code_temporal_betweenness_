#include <string>

template <class V> class FibonacciHeap;

template <class V> struct node {
private:
	node<V>* prev;
	node<V>* next;
	node<V>* child;
	node<V>* parent;
	V value;
	int degree;
	bool marked;
public:
	friend class FibonacciHeap<V>;
	node<V>* getPrev() {return prev;}
	node<V>* getNext() {return next;}
	node<V>* getChild() {return child;}
	node<V>* getParent() {return parent;}
	V getValue() {return value;}
	bool isMarked() {return marked;}

	bool hasChildren() {return child;}
	bool hasParent() {return parent;}
};

template <class V> class FibonacciHeap {
protected:
	node<V>* heap;
public:

	FibonacciHeap() {
		heap=_empty();
	}
	virtual ~FibonacciHeap() {
		if(heap) {
			_deleteAll(heap);
		}
	}
	node<V>* insert(V value) {
		node<V>* ret=_singleton(value);
		heap=_merge(heap,ret);
		return ret;
	}
	void merge(FibonacciHeap& other) {
		heap=_merge(heap,other.heap);
		other.heap=_empty();
	}

	bool isEmpty() {
		return heap==NULL;
	}

	V getMinimum() {
		return heap->value;
	}

	V removeMinimum() {
		node<V>* old=heap;
		heap=_removeMinimum(heap);
		V ret=old->value;
		delete old;
		return ret;
	}

	void decreaseKey(node<V>* n,V value) {
		heap=_decreaseKey(heap,n,value);
	}

	node<V>* find(V value) {
		return _find(heap,value);
	}

	void display() const {	// function code adapted from GO code just below C++
		node<V>* p = heap;
		if (p == NULL) {
			std::cout << "The Heap is Empty" << std::endl;
			return;
		}
		std::cout << "The root nodes of Heap are: " << std::endl;
		_display_tree(heap, "");
		std::cout << std::endl;
	}

private:
	node<V>* _empty() {
		return NULL;
	}

	node<V>* _singleton(V value) {
		node<V>* n=new node<V>;
		n->value=value;
		n->prev=n->next=n;
		n->degree=0;
		n->marked=false;
		n->child=NULL;
		n->parent=NULL;
		return n;
	}

	node<V>* _merge(node<V>* a,node<V>* b) {
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

	void _deleteAll(node<V>* n) {
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
	
	void _addChild(node<V>* parent,node<V>* child) {
		child->prev=child->next=child;
		child->parent=parent;
		parent->degree++;
		parent->child=_merge(parent->child,child);
	}

	void _unMarkAndUnParentAll(node<V>* n) {
		if(n==NULL)return;
		node<V>* c=n;
		do {
			c->marked=false;
			c->parent=NULL;
			c=c->next;
		}while(c!=n);
	}

	node<V>* _removeMinimum(node<V>* n) {
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

	node<V>* _cut(node<V>* heap,node<V>* n) {
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

	node<V>* _decreaseKey(node<V>* heap,node<V>* n,V value) {
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

	node<V>* _find(node<V>* heap,V value) {
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

	void _display_tree(node<V>* n, std::string pre) const {
		std::string pc = "│  ";
		node<V>* x = n;
		do {
			if (x->next != n) {
				std::cout << pre << "├─";
			} else {
				std::cout << pre << "└─";
				pc = "   ";
			}
			if (x->child == nullptr) {
				std::cout << "─" << x->value << std::endl;
			} else {
				std::cout << "┐" << x->value << std::endl;
				_display_tree(x->child, pre + pc);
			}
			//		std::cout << std::endl;
			x = x->next;
		} while (x != n);
	}

};

// /*
//  * main() for testing constructor, getMinimum(), display(), removeMinimum(), decreaseKey(), isEmpty()
//  */
// int main(int argc, char** argv) {

// 	FibonacciHeap<int> fh;

// 	fh.insert(23);
// 	fh.insert(7);
// 	fh.insert(21);
// 	fh.insert(3);
// 	fh.insert(17);
// 	fh.insert(24);
// 	fh.insert(18);
// 	fh.insert(52);
// 	fh.insert(38);
// 	fh.insert(30);
// 	fh.insert(26);
// 	fh.insert(46);
// 	node<int>* n = fh.insert(39);
// 	node<int>* m = fh.insert(41);
// 	fh.insert(35);

// 	cout << "Heap Minimum: " << fh.getMinimum() << endl;
// 	cout << "The Heap is: " << endl;

// 	fh.display();
// 	cout << "Heap Minimum Extracted: " << fh.removeMinimum() << endl;
// 	fh.display();

// 	cout << "de: " << n->getValue() << " para: 5" << endl;
// 	fh.decreaseKey(n, 5);

// 	cout << "Heap Minimum: " << fh.getMinimum() << endl;
// 	fh.display();

// 	cout << "de: " << m->getValue() << " para: 2" << endl;
// 	fh.decreaseKey(m, 2);

// 	while (!fh.isEmpty()) {
// 		cout << "Heap Minimum Extracted: " << fh.removeMinimum() << endl;
// 		fh.display();
// 	}

// 	return 0;
// }

// // Operations on a Fibonacci heap in C++

// #include <cmath>
// #include <cstdlib>
// #include <iostream>

// using namespace std;

// // Node creation
// struct node {
//   int n;
//   int degree;
//   node *parent;
//   node *child;
//   node *left;
//   node *right;
//   char mark;

//   char C;
// };

// // Implementation of Fibonacci heap
// class FibonacciHeap {
//    private:
//   int nH;

//   node *H;

//    public:
//   node *InitializeHeap();
//   void Fibonnaci_link(node *, node *, node *);
//   node *Create_node(int);
//   node *Insert(node *, node *);
//   node *Union(node *, node *);
//   node *Extract_Min(node *);
//   int Consolidate(node *);
//   int Display(node *);
//   node *Find(node *, int);
//   int Decrease_key(node *, int, int);
//   int Delete_key(node *, int);
//   void Cut(node *, node *, node *);
//   void Cascase_cut(node *, node *);
//   FibonacciHeap() { H = InitializeHeap(); }
// };

// // Initialize heap
// node *FibonacciHeap::InitializeHeap() {
//   node *np;
//   np = NULL;
//   return np;
// }

// // Create node
// node *FibonacciHeap::Create_node(int value) {
//   node *x = new node;
//   x->n = value;
//   return x;
// }

// // Insert node
// node *FibonacciHeap::Insert(node *H, node *x) {
//   x->degree = 0;
//   x->parent = NULL;
//   x->child = NULL;
//   x->left = x;
//   x->right = x;
//   x->mark = 'F';
//   x->C = 'N';
//   if (H != NULL) {
//     (H->left)->right = x;
//     x->right = H;
//     x->left = H->left;
//     H->left = x;
//     if (x->n < H->n)
//       H = x;
//   } else {
//     H = x;
//   }
//   nH = nH + 1;
//   return H;
// }

// // Create linking
// void FibonacciHeap::Fibonnaci_link(node *H1, node *y, node *z) {
//   (y->left)->right = y->right;
//   (y->right)->left = y->left;
//   if (z->right == z)
//     H1 = z;
//   y->left = y;
//   y->right = y;
//   y->parent = z;

//   if (z->child == NULL)
//     z->child = y;

//   y->right = z->child;
//   y->left = (z->child)->left;
//   ((z->child)->left)->right = y;
//   (z->child)->left = y;

//   if (y->n < (z->child)->n)
//     z->child = y;
//   z->degree++;
// }

// // Union Operation
// node *FibonacciHeap::Union(node *H1, node *H2) {
//   node *np;
//   node *H = InitializeHeap();
//   H = H1;
//   (H->left)->right = H2;
//   (H2->left)->right = H;
//   np = H->left;
//   H->left = H2->left;
//   H2->left = np;
//   return H;
// }

// // Display the heap
// int FibonacciHeap::Display(node *H) {
//   node *p = H;
//   if (p == NULL) {
//     cout << "Empty Heap" << endl;
//     return 0;
//   }
//   cout << "Root Nodes: " << endl;

//   do {
//     cout << p->n;
//     p = p->right;
//     if (p != H) {
//       cout << "-->";
//     }
//   } while (p != H && p->right != NULL);
//   cout << endl;
// }

// // Extract min
// node *FibonacciHeap::Extract_Min(node *H1) {
//   node *p;
//   node *ptr;
//   node *z = H1;
//   p = z;
//   ptr = z;
//   if (z == NULL)
//     return z;

//   node *x;
//   node *np;

//   x = NULL;

//   if (z->child != NULL)
//     x = z->child;

//   if (x != NULL) {
//     ptr = x;
//     do {
//       np = x->right;
//       (H1->left)->right = x;
//       x->right = H1;
//       x->left = H1->left;
//       H1->left = x;
//       if (x->n < H1->n)
//         H1 = x;

//       x->parent = NULL;
//       x = np;
//     } while (np != ptr);
//   }

//   (z->left)->right = z->right;
//   (z->right)->left = z->left;
//   H1 = z->right;

//   if (z == z->right && z->child == NULL)
//     H = NULL;

//   else {
//     H1 = z->right;
//     Consolidate(H1);
//   }
//   nH = nH - 1;
//   return p;
// }

// // Consolidation Function
// int FibonacciHeap::Consolidate(node *H1) {
//   int d, i;
//   float f = (log(nH)) / (log(2));
//   int D = f;
//   node** A = new node*[D];
//   //node *A[D];

//   for (i = 0; i <= D; i++)
//     A[i] = NULL;

//   node *x = H1;
//   node *y;
//   node *np;
//   node *pt = x;

//   do {
//     pt = pt->right;

//     d = x->degree;

//     while (A[d] != NULL)

//     {
//       y = A[d];

//       if (x->n > y->n)

//       {
//         np = x;

//         x = y;

//         y = np;
//       }

//       if (y == H1)
//         H1 = x;
//       Fibonnaci_link(H1, y, x);
//       if (x->right == x)
//         H1 = x;
//       A[d] = NULL;
//       d = d + 1;
//     }

//     A[d] = x;
//     x = x->right;

//   }

//   while (x != H1);
//   H = NULL;
//   for (int j = 0; j <= D; j++) {
//     if (A[j] != NULL) {
//       A[j]->left = A[j];
//       A[j]->right = A[j];
//       if (H != NULL) {
//         (H->left)->right = A[j];
//         A[j]->right = H;
//         A[j]->left = H->left;
//         H->left = A[j];
//         if (A[j]->n < H->n)
//           H = A[j];
//       } else {
//         H = A[j];
//       }
//       if (H == NULL)
//         H = A[j];
//       else if (A[j]->n < H->n)
//         H = A[j];
//     }
//   }
// }

// // Decrease Key Operation
// int FibonacciHeap::Decrease_key(node *H1, int x, int k) {
//   node *y;
//   if (H1 == NULL) {
//     cout << "The Heap is Empty" << endl;
//     return 0;
//   }
//   node *ptr = Find(H1, x);
//   if (ptr == NULL) {
//     cout << "Node not found in the Heap" << endl;
//     return 1;
//   }

//   if (ptr->n < k) {
//     cout << "Entered key greater than current key" << endl;
//     return 0;
//   }
//   ptr->n = k;
//   y = ptr->parent;
//   if (y != NULL && ptr->n < y->n) {
//     Cut(H1, ptr, y);
//     Cascase_cut(H1, y);
//   }

//   if (ptr->n < H->n)
//     H = ptr;

//   return 0;
// }

// // Cutting Function
// void FibonacciHeap::Cut(node *H1, node *x, node *y)

// {
//   if (x == x->right)
//     y->child = NULL;
//   (x->left)->right = x->right;
//   (x->right)->left = x->left;
//   if (x == y->child)
//     y->child = x->right;
//   y->degree = y->degree - 1;
//   x->right = x;
//   x->left = x;
//   (H1->left)->right = x;
//   x->right = H1;
//   x->left = H1->left;
//   H1->left = x;
//   x->parent = NULL;
//   x->mark = 'F';
// }

// // Cascade cut
// void FibonacciHeap::Cascase_cut(node *H1, node *y) {
//   node *z = y->parent;
//   if (z != NULL) {
//     if (y->mark == 'F') {
//       y->mark = 'T';
//     } else

//     {
//       Cut(H1, y, z);
//       Cascase_cut(H1, z);
//     }
//   }
// }

// // Search function
// node *FibonacciHeap::Find(node *H, int k) {
//   node *x = H;
//   x->C = 'Y';
//   node *p = NULL;
//   if (x->n == k) {
//     p = x;
//     x->C = 'N';
//     return p;
//   }

//   if (p == NULL) {
//     if (x->child != NULL)
//       p = Find(x->child, k);
//     if ((x->right)->C != 'Y')
//       p = Find(x->right, k);
//   }

//   x->C = 'N';
//   return p;
// }

// // Deleting key
// int FibonacciHeap::Delete_key(node *H1, int k) {
//   node *np = NULL;
//   int t;
//   t = Decrease_key(H1, k, -5000);
//   if (!t)
//     np = Extract_Min(H);
//   if (np != NULL)
//     cout << "Key Deleted" << endl;
//   else
//     cout << "Key not Deleted" << endl;
//   return 0;
// }

// int main() {
//   int n, m, l;
//   FibonacciHeap fh;
//   node *p;
//   node *H;
//   H = fh.InitializeHeap();

//   p = fh.Create_node(7);
//   H = fh.Insert(H, p);
//   p = fh.Create_node(3);
//   H = fh.Insert(H, p);
//   p = fh.Create_node(17);
//   H = fh.Insert(H, p);
//   p = fh.Create_node(24);
//   H = fh.Insert(H, p);

//   fh.Display(H);

//   p = fh.Extract_Min(H);
//   if (p != NULL)
//     cout << "The node with minimum key: " << p->n << endl;
//   else
//     cout << "Heap is empty" << endl;

//   m = 26;
//   l = 16;
//   fh.Decrease_key(H, m, l);

//   m = 16;
//   fh.Delete_key(H, m);
// }
