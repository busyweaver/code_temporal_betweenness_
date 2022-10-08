// A C++ program to demonstrate common Binary Heap Operations
#include<iostream>
#include<climits>
#include<map>
#include<numeric>
#include <limits>
using namespace std;
std::ostream& operator<< (std::ostream& os, const std::pair<double,std::pair<int,int>> & lhs)
{
  os << "("<<lhs.first
     << ", "
     << lhs.second.first
     << ", "
     << lhs.second.second << ")";
  return os; 
}


bool less_pair( std::pair<double,std::pair<int,int>>  const a, std::pair<double,std::pair<int,int>>  const b)
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
bool great_pair( std::pair<double,std::pair<int,int>>  const a, std::pair<double,std::pair<int,int>>  const b)
{
  return ((!(less_pair(a,b)) )  && (!(equal_pair(a,b))));
}
// Prototype of a utility function to swap two integers

// A class for Min Heap
class MinHeap
{

	int capacity; // maximum possible size of min heap


public:
  std::pair<double,std::pair<int,int>> *harr; // pointer to array of elements in heap
  int heap_size; // Current number of elements in min heap
  std::map<std::pair<double,std::pair<int,int>>, int> index_elem;
	// Constructor
	MinHeap(int capacity);

	// to heapify a subtree with the root at given index
	void MinHeapify(int );

	int parent(int i) { return (i-1)/2; }

	// to get index of left child of node at index i
	int left(int i) { return (2*i + 1); }

	// to get index of right child of node at index i
	int right(int i) { return (2*i + 2); }

	// to extract the root which is the minimum element
	std::pair<double,std::pair<int,int>> extractMin();

	// Decreases key value of key at index i to new_val
	void decreaseKey(int i, std::pair<double,std::pair<int,int>> new_val);

	// Returns the minimum key (key at root) from min heap
	std::pair<double,std::pair<int,int>> getMin() { return harr[0]; }

	// Deletes a key stored at index i
	void deleteKey(int i);

	// Inserts a new key 'k'
	void insertKey(std::pair<double,std::pair<int,int>> k);
  void swap(std::pair<double,std::pair<int,int>> *x, std::pair<double,std::pair<int,int>> *y);
  void display();
};

// Constructor: Builds a heap from a given array a[] of given size
MinHeap::MinHeap(int cap)
{
	heap_size = 0;
	capacity = cap;
	harr = new std::pair<double,std::pair<int,int>>[cap];

}

// Inserts a new key 'k'
void MinHeap::insertKey(std::pair<double,std::pair<int,int>> k)
{
	if (heap_size == capacity)
	{
		cout << "\nOverflow: Could not insertKey\n";
		return;
	}

	// First insert the new key at the end
	heap_size++;
	int i = heap_size - 1;
	harr[i] = k;
  index_elem[k] = i;

	// Fix the min heap property if it is violated
	while (i != 0 && great_pair(harr[parent(i)] , harr[i]) )
	{
	swap(&harr[i], &harr[parent(i)]);
	i = parent(i);
	}
}

// Decreases value of key at index 'i' to new_val. It is assumed that
// new_val is smaller than harr[i].
void MinHeap::decreaseKey(int i, std::pair<double,std::pair<int,int>> new_val)
{
  index_elem[new_val] = index_elem[harr[i]];
  index_elem.erase(harr[i]);
	harr[i] = new_val;

	while (i != 0 && great_pair(harr[parent(i)] , harr[i]))
	{
	swap(&harr[i], &harr[parent(i)]);
	i = parent(i);
	}
}

// Method to remove minimum element (or root) from min heap
std::pair<double,std::pair<int,int>> MinHeap::extractMin()
{
	if (heap_size <= 0)
    {
      std::pair<double,std::pair<int,int>> p;
      p.first = std::numeric_limits<double>::infinity();
      p.second.first = -1;
      p.second.second = -1;
      return p;
    }

	if (heap_size == 1)
	{
		heap_size--;
    index_elem.erase(harr[0]);
		return harr[0];

	}

	// Store the minimum value, and remove it from heap
  index_elem.erase(harr[0]);
	std::pair<double,std::pair<int,int>> root = harr[0];
	harr[0] = harr[heap_size-1];
  index_elem[harr[0]] = 0;
	heap_size--;
	MinHeapify(0);

	return root;
}


// This function deletes key at index i. It first reduced value to minus
// infinite, then calls extractMin()
void MinHeap::deleteKey(int i)
{
  std::pair<double,std::pair<int,int>> p;
  p.first = std::numeric_limits<double>::infinity();
  p.second.first = -1;
  p.second.second = -1;
	decreaseKey(i, p);
	extractMin();
}

// A recursive method to heapify a subtree with the root at given index
// This method assumes that the subtrees are already heapified
void MinHeap::MinHeapify(int i)
{
	int l = left(i);
	int r = right(i);
	int smallest = i;
	if (l < heap_size && less_pair(harr[l] , harr[i]))
		smallest = l;
	if (r < heap_size && less_pair(harr[r] , harr[smallest]))
		smallest = r;
	if (smallest != i)
	{
		swap(&harr[i], &harr[smallest]);
		MinHeapify(smallest);
	}

}

void MinHeap::swap(std::pair<double,std::pair<int,int>> *x, std::pair<double,std::pair<int,int>> *y)
{
  std::pair<double,std::pair<int,int>> temp = *x;
  int tmp;
  tmp = index_elem[*x];
  index_elem[*x] = index_elem[*y];
  index_elem[*y] = tmp;
  *x = *y;
  *y = temp;
}

void MinHeap::display()
{
  printf("\nheap\n");
  for(int i = 0; i < heap_size ; i++)
    std::cout << harr[i] << " ";
  printf("\nend heap\n");
}

// A utility function to swap two elements


// Driver program to test above functions
/* int main() */
/* { */
/* 	MinHeap h(11); */
/* 	h.insertKey(3); */
/* 	h.insertKey(2); */
/* 	h.deleteKey(1); */
/* 	h.insertKey(15); */
/* 	h.insertKey(5); */
/* 	h.insertKey(4); */
/* 	h.insertKey(45); */
/* 	cout << h.extractMin() << " "; */
/* 	cout << h.getMin() << " "; */
/* 	h.decreaseKey(2, 1); */
/* 	cout << h.getMin(); */
/* 	return 0; */
/* } */
