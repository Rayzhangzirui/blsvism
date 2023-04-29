#ifndef HEAP_H
#define HEAP_H

#include "util.h"

struct HeapElt
{
   int index[3];
   HeapElt* child[2];
   HeapElt* parent;
   int num;
};

struct HeapStruct
{
   HeapElt *head;
   HeapElt *tail;
   HeapElt ****loc; // 3d grid to HeapElt*
   double ***value;
   char ***tube;
   static const int dim = 3;
};

HeapElt ****heapmatrix(int row, int col, int fr);

// free heap matrix, HeapElt 3d
void free_heapmatrix(HeapElt ****m, int row, int col, int fr);


void addtoheap(HeapStruct &heap, double value, int *index);
HeapElt* getheapelt(HeapStruct heap, int num);
void fixheapeltup(HeapStruct &heap, HeapElt *fix);
HeapElt* fixheapeltempty(HeapStruct &heap, HeapElt *fix);
void fixheapeltdelete(HeapStruct &heap, HeapElt *del);
void fixheapeltreplace(HeapStruct &heap, double value, int *index);
char betterheapval(double val1, double val2);
void heaptotals(HeapStruct heap);

#endif



