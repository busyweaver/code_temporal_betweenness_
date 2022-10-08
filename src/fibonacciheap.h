/*Copyright (c) 2010, Robin Message <Robin.Message@cl.cam.ac.uk>
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Univsersity of Cambridge nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF CAMBRIDGE OR ROBIN MESSAGE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <iostream>     // std::cout
#include <string>
std::ostream& operator<< (std::ostream& os, const std::pair<double,std::pair<int,int>> & lhs)
{
  os << lhs.first
     << ", "
     << lhs.second.first
     << ", "
     << lhs.second.second;
  return os; 
}


bool cmp_pair( std::pair<double,std::pair<int,int>>  const a, std::pair<double,std::pair<int,int>>  const b)
{
  if (a.first < b.first)
    {
      return true;
    }
  else if (a.first == b.first)
    {
      if (a.second.first < b.second.first)
        return true;
      else if(a.second.first == b.second.first)
        return (a.second.second < b.second.second);
      else
        return false;
    }
  else
    return false;
}
bool equal_pair( std::pair<double,std::pair<int,int>>  const a, std::pair<double,std::pair<int,int>>  const b)
{
  if (a.first == b.first && a.second.first == b.second.first && a.second.second == b.second.second)
    return true;
  return false;
}
class FibonacciHeap;

 struct node {
private:
	node* prev;
	node* next;
	node* child;
	node* parent;
	std::pair<double,std::pair<int,int>> value;
	int degree;
	bool marked;
public:
	friend class FibonacciHeap;
	node* getPrev() {return prev;}
	node* getNext() {return next;}
	node* getChild() {return child;}
	node* getParent() {return parent;}
   std::pair<double,std::pair<int,int>> getValue() {return value;}
	bool isMarked() {return marked;}

	bool hasChildren() {return child;}
	bool hasParent() {return parent;}
};

 class FibonacciHeap {
protected:
	node* heap;
public:

	FibonacciHeap() {
		heap=_empty();
	}
	virtual ~FibonacciHeap() {
		if(heap) {
			_deleteAll(heap);
		}
	}
	node* insert(std::pair<double,std::pair<int,int>> value) {
		node* ret=_singleton(value);
		heap=_merge(heap,ret);
		return ret;
	}
	void merge(FibonacciHeap& other) {
		heap=_merge(heap,other.heap);
		other.heap=_empty();
	}

	bool isEmpty() {
		return heap==nullptr;
	}

	std::pair<double,std::pair<int,int>> getMinimum() {
		return heap->value;
	}

	std::pair<double,std::pair<int,int>> removeMinimum() {
    printf("remove\n");
		node* old=heap;
		heap=_removeMinimum(heap);
		//std::pair<double,std::pair<int,int>> ret=old->value;
    std::pair<double,std::pair<int,int>> ret;
    ret.first = old->value.first;
    ret.second.first = old->value.second.first;
    ret.second.second = old->value.second.second;
    std::cout << ret << "\n";
		delete old;
		return ret;
	}

	void decreaseKey(node* n,std::pair<double,std::pair<int,int>> value) {
		heap=_decreaseKey(heap,n,value);
	}

	node* find(std::pair<double,std::pair<int,int>> value) {
		return _find(heap,value);
	}

	void display() const {	// function code adapted from GO code just below C++
		node* p = heap;
		if (p == nullptr) {
			std::cout << "The Heap is Empty" << std::endl;
			return;
		}
		std::cout << "The root nodes of Heap are: " << std::endl;
		_display_tree(heap, "");
		std::cout <<  " "   << std::endl;
	}

private:
	node* _empty() {
		return nullptr;
	}

	node* _singleton(std::pair<double,std::pair<int,int>> value) {
		node* n=new node;
		n->value.first=value.first;
    n->value.second.first=value.second.first;
    n->value.second.second=value.second.second;
		n->prev=n->next=n;
		n->degree=0;
		n->marked=false;
		n->child=nullptr;
		n->parent=nullptr;
		return n;
	}

	node* _merge(node* a,node* b) {
		if(a==nullptr)return b;
		if(b==nullptr)return a;
		if( (!(cmp_pair(a->value,b->value))) && (!(equal_pair(a->value,b->value)))   ) {
			node* temp=a;
			a=b;
			b=temp;
		}
		node* an=a->next;
		node* bp=b->prev;
		a->next=b;
		b->prev=a;
		an->prev=bp;
		bp->next=an;
		return a;
	}

	void _deleteAll(node* n) {
		if(n!=nullptr) {
			node* c=n;
			do {
				node* d=c;
				c=c->next;
				_deleteAll(d->child);
				delete d;
			} while(c!=n);
		}
	}
	
	void _addChild(node* parent,node* child) {
		child->prev=child->next=child;
		child->parent=parent;
		parent->degree++;
		parent->child=_merge(parent->child,child);
	}

	void _unMarkAndUnParentAll(node* n) {
		if(n==nullptr)return;
		node* c=n;
		do {
			c->marked=false;
			c->parent=nullptr;
			c=c->next;
		}while(c!=n);
	}

	node* _removeMinimum(node* n) {
    std::cout << "_remve \n";
		_unMarkAndUnParentAll(n->child);
		if(n->next==n) {
			n=n->child;
		} else {
			n->next->prev=n->prev;
			n->prev->next=n->next;
			n=_merge(n->next,n->child);
		}
		if(n==nullptr)return n;
		node* trees[64]={nullptr};
		
		while(true) {
			if(trees[n->degree]!=nullptr) {
				node* t=trees[n->degree];
				if(t==n)break;
				trees[n->degree]=nullptr;
				t->prev->next=t->next;
				t->next->prev=t->prev;
				if(cmp_pair(n->value,t->value)) {
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
		node* min=n;
		do {
      std::cout << n->value << " deuxieme "<< min->value << " cmp "<<cmp_pair(n->value,min->value)<< "\n";
			if(cmp_pair(n->value,min->value))min=n;
			n=n->next;
		} while(n!=n);
		return min;
	}

	node* _cut(node* heap,node* n) {
		if(n->next==n) {
			n->parent->child=nullptr;
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

	node* _decreaseKey(node* heap,node* n,std::pair<double,std::pair<int,int>> value) {
		if(cmp_pair(n->value,value))return heap;
		n->value=value;
		node* parent = n->parent;
		if(parent != nullptr && cmp_pair(n->value , parent->value)) {
			heap=_cut(heap,n);
			n->parent=nullptr;
			while(parent!=nullptr && parent->marked) {
				heap=_cut(heap,parent);
				n=parent;
				parent=n->parent;
				n->parent=nullptr;
			}
			if(parent!=nullptr && parent->parent!=nullptr)parent->marked=true;
			if (cmp_pair(n->value , heap->value))heap = n;
		}
		return heap;
	}

	node* _find(node* heap,std::pair<double,std::pair<int,int>> value) {
		node* n=heap;
		if(n==nullptr)return nullptr;
		do {
			if(equal_pair(n->value,value))return n;
			node* ret=_find(n->child,value);
			if(ret)return ret;
			n=n->next;
		}while(n!=heap);
		return nullptr;
	}

	void _display_tree(node* n, std::string pre) const {
		std::string pc = "│  ";
		node* x = n;
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

