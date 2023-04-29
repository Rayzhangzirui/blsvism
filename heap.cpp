#include "heap.h"

HeapElt ****heapmatrix(int row, int col, int fr)
{
   int i, j, k;
   HeapElt ****m;
   int err = 0;

   m = new HeapElt***[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new HeapElt**[col+1];
      if (!m[i] && !err)
       err = 1;
      for (j = 0; j <= col; j++)
      {
       m[i][j] = new HeapElt*[fr+1];
       if (!m[i] && !err)
         err = 1;
       for (k = 0; k <= fr; k++)
         m[i][j][k] = NULL;
      }
   }
   return m;
}

// free heap matrix, HeapElt 3d
void free_heapmatrix(HeapElt ****m, int row, int col, int fr)
{
   int i, j, k;

   for (i = 0; i <= row; i++){
      for (j = 0; j <= col; j++){
         for (k = 0; k <= fr; k++){
            delete [] (m[i][j][k]); 
         }
         delete [] (m[i][j]);
      }
      delete [] (m[i]);
   }
   delete [] m;
}

/*
add element at the end of the heap, and then bubble up
store value in heap.value
add pointer from heap.loc to element
In the whole algorithm, the structure of the heapElt does not change,
only the values/index in the heapElt is changed
*/
void addtoheap(HeapStruct &heap, double value, int *index)
{
	HeapElt *newelt, *theparent;
	int i, last;

	// last is the current number of elements in the heap, starting from 1
	// last+1 is the next slot
	if (heap.tail != NULL)
		last = (*(heap.tail)).num;
	else
		last = 0;

	newelt = getheapelt(heap,last+1);

	// new element
	if (newelt == NULL)
	{
		if (last > 0)
		{
		 theparent = getheapelt(heap,(int) ((last+1)/2));//theparent points to parent heapElt
		 
		 //add child, left or right, depending of last
		 // if (last+1) is even, left child
		 (*theparent).child[(last+1)%2] = new HeapElt;
		 newelt = (*theparent).child[(last+1)%2];//go to child, connect to parent
		 (*newelt).parent = theparent;
		}
		else // if at here, then we are inserting first element to an empty heap
		{
		 heap.head = new HeapElt;
		 newelt = heap.head;
		 (*newelt).parent = NULL;
		}
		(*newelt).num = last+1;
		(*newelt).child[0] = NULL;
		(*newelt).child[1] = NULL;
	}
	for (i = 0; i < heap.dim; i++)
		(*newelt).index[i] = index[i];
	setvalarray(heap.value,index,value);
	setvalarray(heap.loc,index,newelt);
	heap.tail = newelt;

	fixheapeltup(heap,newelt);
}


// navigate to num-th element of heap, head is num 1. let j be the layer num
// at the root, j=1, num=1.
HeapElt* getheapelt(HeapStruct heap, int num)
{
	HeapElt *current;
	int n, k;

	if (num <= 0)
		return NULL;

	for (n = 1; n <= num; n*= 2);// keep multiplying by 2, untial n = 2^(j+1), first element in the next layer below num
	n /= 2; //n = 2^j, first element same layer as num

	current = heap.head;
	k = num;
	// the loop tell you to go left or go right starting from the root to reach num
	// at the j row, there are 2^j element, by compare distance between num and the first element in the same row
	// if less than 2^j/2, goes left, otherwise, goes right
	while (n > 1 && current != NULL)
	{
		if (k >= n)
		 k -= n;
		n /= 2;
		if (k >= n)
		 current = (*current).child[1];
		else
		 current = (*current).child[0];
	}

	return current;
}

/*
sift up: compare value of child and parent
swap index of child and parent, also swap heap.loc
*/ 
void fixheapeltup(HeapStruct &heap, HeapElt *fix)
{
	HeapElt *parent, *child, temp;
	char change = 1;
	int i;

	child = fix;
	parent = (*child).parent;
	while (change && parent != NULL)
	{
		if (betterheapval(evalarray(heap.value,(*child).index),
						evalarray(heap.value,(*parent).index)))//if child has smaller value then parent, swap
		{
		 for (i = 0; i < heap.dim; i++) 
			temp.index[i] = (*parent).index[i];
		 for (i = 0; i < heap.dim; i++) 
			(*parent).index[i] = (*child).index[i];
		 setvalarray(heap.loc,(*parent).index,parent);
		 for (i = 0; i < heap.dim; i++) 
			(*child).index[i] = temp.index[i];
		 setvalarray(heap.loc,(*child).index,child);
		 child = parent;
		 parent = (*parent).parent;
		}
		else
		 change = 0;
	}
}

// remove HeapElt* fix, sift up child with smaller value, should delete fix at the end
// fix will go the the leave of the tree
// empty = disconnect heap.loc to HeapElt
HeapElt* fixheapeltempty(HeapStruct &heap, HeapElt *fix)
{
	HeapElt *current, *child[2];
	int i, r;

	setvalarray<HeapElt*>(heap.loc,(*fix).index, NULL);

	current = fix;
	child[0] = (*current).child[0];
	child[1] = (*current).child[1];
	while (child[0] != NULL && (*(child[0])).index[0] >= 0)
	{
		// no right child OR right child is invalid OR right child has worse value (larger du)
		if (child[1] == NULL || (*(child[1])).index[0] < 0 ||
			betterheapval(evalarray(heap.value,(*(child[0])).index),
						evalarray(heap.value,(*(child[1])).index)))
		 r = 0; //child element with better value
		else
		 r = 1;
		for (i = 0; i < heap.dim; i++)
		 (*current).index[i] = (*(child[r])).index[i];//copy child[r] to current
		setvalarray<HeapElt*>(heap.loc,(*current).index,current);//update heap.loc point to current
		current = child[r];//point current to child[r]
		child[0] = (*current).child[0];
		child[1] = (*current).child[1];
	}

	return current;
}

/*
delete element from the heap, when energy positive, i.e. when flip the top
*/
void fixheapeltdelete(HeapStruct &heap, HeapElt *del)
{
	HeapElt *current;
	int i;
	HeapElt *temp = heap.tail;

	//empty del, bubble value up
	current = fixheapeltempty(heap,del);
	if (current != heap.tail)
	{
		// if current lands in some leaf that is not the tail
		// swap index of current and tail, and sift current up
		// set index to be (-1,-1,-1)
		// and  move tail pointer before by one.
		for (i = 0; i < heap.dim; i++)   
		 (*current).index[i] = (*(heap.tail)).index[i];
		setvalarray(heap.loc,(*current).index,current);//copy tail to current
		for (i = 0; i < heap.dim; i++)   
		 (*(heap.tail)).index[i] = -1; //set tail to have invalid index
		
		heap.tail = getheapelt(heap,(*(heap.tail)).num-1);

		fixheapeltup(heap,current);
	}
	else
	{
		// set tail index to be (-1,-1,-1)
		// and  move tail pointer before by one.
		for (i = 0; i < heap.dim; i++)   
		 (*(heap.tail)).index[i] = -1;
		heap.tail = getheapelt(heap,(*(heap.tail)).num-1);
	}
	
	//delete tail element. The tail element cannot be accessed from heap.loc anymore
	if (temp->parent){
		if ((temp->parent)->child[0] == temp){
		temp->parent->child[0] = NULL;
		}else{
			temp->parent->child[1] = NULL;
		}	
		delete temp;
	}
	
	

}

// called when modify du of elements within kernel rad
// change value of heapElt: modify heap.value,
// empty heapElt(index), so it drop down to the leaf
// then reconnect heap.value to heapElt(index), sift up
void fixheapeltreplace(HeapStruct &heap, double value, int *index)
{
	HeapElt *replace, *current, temp;
	int i;

	replace = evalarray(heap.loc,index);

	for (i = 0; i < heap.dim; i++)   
		temp.index[i] = (*replace).index[i];
	setvalarray(heap.value,temp.index,value); //temp copy value and index

	current = fixheapeltempty(heap,replace);//remove replace, bubble child up, current reaches bottom, 
	for (i = 0; i < heap.dim; i++)   
		(*current).index[i] = temp.index[i];
	setvalarray(heap.loc,(*current).index,current);// put temp at current
	fixheapeltup(heap,current);// fix heap
}

char betterheapval(double val1, double val2)
{
	if (val1 < val2)
		return 1;
	else
		return 0;
}


void heaptotals(HeapStruct heap)
{
	int r, num;
	HeapElt *current = heap.tail;
	if (current != NULL)
	{
		num = (*current).num;
		current = getheapelt(heap,num+1);
	}
	else
	{
		num = 0;
		current = heap.head;
	}

	for (r = num+1; current != NULL; 
		r++, current = getheapelt(heap,r));
	
	cout << "Heap has " << num << " active nodes and " << r-1 << " total nodes" 
		<< endl;
}