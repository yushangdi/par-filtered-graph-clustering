#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/primitives.h"

#include "symmatrix.h"


namespace IO {
  using namespace std;

  auto is_space = [] (char c) {
    switch (c)  {
    case '\r': 
    case '\t': 
    case '\n': 
    case 0:
    case ' ' : return true;
    default : return false;
    }
  };


  // modified from pbbsbench  https://github.com/cmuparlay/pbbsbench/blob/37df3e14c7f3d738500e06840874a1505944598d/common/geometryIO.h
  parlay::sequence<char> readStringFromFile(char const *fileName) {
    ifstream file (fileName, ios::in | ios::binary | ios::ate);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      abort();
    }
    long end = file.tellg();
    file.seekg (0, ios::beg);
    long n = end - file.tellg();
    parlay::sequence<char> bytes(n, (char) 0);
    file.read (bytes.begin(), n);
    file.close();
    return bytes;
  }
  
  // parallel code for converting a string to word pointers
  // side effects string by setting to null after each word
  template <class Seq>
    parlay::sequence<char*> stringToWords(Seq &Str) {
    size_t n = Str.size();
    
    parlay::parallel_for(0, n, [&] (long i) {
	if (is_space(Str[i])) Str[i] = 0;}); 

    // mark start of words
    auto FL = parlay::tabulate(n, [&] (long i) -> bool {
	return (i==0) ? Str[0] : Str[i] && !Str[i-1];});
    
    // offset for each start of word
    auto Offsets = parlay::pack_index<long>(FL);

    // pointer to each start of word
    auto SA = parlay::tabulate(Offsets.size(), [&] (long j) -> char* {
	return Str.begin() + Offsets[j];});
    
    return SA;
  }

  template <class T, class Seq>
  SymM<T> parseSymMatrix(Seq W, std::size_t n) {
    SymM<T> matrix = SymM<T>(n);
    parlay::parallel_for(0, n, [&](size_t i){
        parlay::parallel_for(i+1, n,[&] (size_t j){
	        matrix.update(i, j, (T)atof(W[i*n + j]));
        });
    });
    return matrix;
  }

  // read a symmatric matrix from file
  template <class T>
  SymM<T> readSymMatrixFromFile(char const *fname, std::size_t n) {
    parlay::sequence<char> S = readStringFromFile(fname);
    parlay::sequence<char*> W = stringToWords(S);
    if (W.size() == 0) {
      cout << "readPointsFromFile empty file" << endl;
      abort();
    }
    if (W.size() % n != 0) {
      cout << "readPointsFromFile wrong file type or wrong dimension" << endl;
      abort();
    }
    return parseSymMatrix<T>(W.cut(0,W.size()), n);
  }

}  // namespace IO