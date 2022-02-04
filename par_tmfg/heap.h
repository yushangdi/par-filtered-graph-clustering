#pragma once

// from fast cluster 
//https://github.com/yushangdi/fastcluster/blob/5f9e688111ca557997b5c24c6fffc2c48b6ea408/src/fastcluster.cpp


#include "parlay/sequence.h"

template<class T>
class binary_min_heap {
  /*
  Class for a binary min-heap. The data resides in an array A. The elements of
  A are not changed but two lists I and R of indices are generated which point
  to elements of A and backwards.
  The heap tree structure is
     H[2*i+1]     H[2*i+2]
         \            /
          \          /
           ≤        ≤
            \      /
             \    /
              H[i]
  where the children must be less or equal than their parent. Thus, H[0]
  contains the minimum. The lists I and R are made such that H[i] = A[I[i]]
  and R[I[i]] = i.
  This implementation is not designed to handle NaN values.
  */

public:
  T* A;
  size_t size;
  size_t* I;
  // size_t* R;

  // no default constructor
  binary_min_heap();

  binary_min_heap(T* A_, const size_t size_, size_t* I_)
    : A(A_), size(size_), I(I_)//, R(R_)
  { // Initialize the lists I and R to the identity. This
    // does not make it a heap. Call heapify afterwards!
    for (size_t i=0; i<size; ++i)
      I[i] = i; //R[i] = 
  }

  ~binary_min_heap() {}

  void heapify() {
    // Arrange the indices I and R so that H[i] := A[I[i]] satisfies the heap
    // condition H[i] < H[2*i+1] and H[i] < H[2*i+2] for each i.
    //
    // Complexity: Θ(size)
    // Reference: Cormen, Leiserson, Rivest, Stein, Introduction to Algorithms,
    // 3rd ed., 2009, Section 6.3 “Building a heap”
    size_t idx;
    for (idx=(size>>1); idx>0; ) {
      --idx;
      update_geq_(idx);
    }
  }

  inline size_t argmin() {
    // Return the minimal element.
    return I[0];
  }

  void heap_pop() {
    // Remove the minimal element from the heap.
#ifdef DEBUG
    if(size==0){
      cout << "size is zero" <<endl;
    }
#endif  
    --size;
    I[0] = I[size];
    // R[I[0]] = 0;
    update_geq_(0);
  }


private:

  void update_geq_ (size_t i) {
    size_t j;
    for ( ; (j=2*i+1)<size; i=j) {
      if ( H(j)>=H(i) ) {
        ++j;
        if ( j>=size || H(j)>=H(i) ) break;
      }
      else if ( j+1<size && H(j+1)<H(j) ) ++j;
      heap_swap(i, j);
    }
  }

  void heap_swap(const size_t i, const size_t j) {
    // Swap two indices.
    size_t tmp = I[i];
    I[i] = I[j];
    I[j] = tmp;
  }

  inline T H(const size_t i) {
    return A[I[i]];
  }

};
