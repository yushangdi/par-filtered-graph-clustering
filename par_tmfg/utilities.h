#pragma once

#include <tuple>
#include <iostream>
#include <fstream> 
#include "parlay/sequence.h"

using namespace std;
// efficient way to read matrix file?
// read a matrix file to sequence A, first two numbers are rows and columns

using vtx=int;
using face=int;
using cliqueT=tuple<vtx,vtx,vtx,vtx>;
using triT=tuple<vtx,vtx,vtx>;

template<class T>
tuple<int, int> readMatrix(std::string filename, parlay::sequence<T> &A){
    // http://www.cplusplus.com/forum/beginner/38349/
    ifstream f(filename); std::size_t m,n;
    f >> m >> n;
    assert(n > 0);
    A = parlay::sequence<T>(m*n);
    for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < n; j++)
    f >> A[i*n + j];
    f.close();
    return make_tuple(m,n);
}

template<class T>
tuple<int, int> readMatrix(std::string filename, parlay::sequence<T> &A, std::size_t n){
    // http://www.cplusplus.com/forum/beginner/38349/
    ifstream f(filename); 
    A = parlay::sequence<T>(n*n);
    for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
    f >> A[i*n + j];
    f.close();
    return make_tuple(n,n);
}

template<class T>
tuple<int, int> readMatrix(std::string filename, parlay::sequence<T> &A, std::size_t n, std::size_t m){
    // http://www.cplusplus.com/forum/beginner/38349/
    ifstream f(filename); 
    A = parlay::sequence<T>(n*m);
    for (std::size_t i = 0; i < n; i++){
    for (std::size_t j = 0; j < m; j++){
    f >> A[i*m + j];
    }}
    f.close();
    return make_tuple(n,m);
}

// struct hash_vertex_to_faces {
//   using eType = tuple<vtx,long,int>; // vertex, index into its array, number of faces
//   using kType = vtx;
//   eType empty() { return make_tuple(-1,-1,0); }
//   kType getKey(eType v) { return get<0>(v); }
//   size_t hash(kType v) { return static_cast<size_t>(hash64(v)); }
//   int cmp(kType v, kType b) { return (v > b) ? 1 : ((v == b) ? 0 : -1); }
//   bool replaceQ(eType, eType) { return 0; }
//   eType update(eType v, eType) { return v; }
//   bool cas(eType* p, eType o, eType n) {
//     // TODO: Make this use atomics properly. This is a quick
//     // fix to get around the fact that the hashtable does
//     // not use atomics. This will break for types that
//     // do not inline perfectly inside a std::atomic (i.e.,
//     // any type that the standard library chooses to lock)
//     return std::atomic_compare_exchange_strong_explicit(
//       reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
//   }
// };
// using tableT=parlay::hashtable<hash_vertex_to_faces>;


//https://stackoverflow.com/questions/6245735/pretty-print-stdtuple
// helper function to print a tuple of any size
template<class Tuple, std::size_t N>
struct TuplePrinter {
    static void print(const Tuple& t) 
    {
        TuplePrinter<Tuple, N-1>::print(t);
        std::cout << ", " << std::get<N-1>(t);
    }
};

template<class Tuple>
struct TuplePrinter<Tuple, 1> {
    static void print(const Tuple& t) 
    {
        std::cout << std::get<0>(t);
    }
};

template<class... Args>
void print(const std::tuple<Args...>& t) 
{
    std::cout << "(";
    TuplePrinter<decltype(t), sizeof...(Args)>::print(t);
    std::cout << ")\n";
}
// end helper function