#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/primitives.h"

// #include "geometry.h"
#include "symmatrix.h"


namespace research_graph {
namespace in_memory {
namespace IO {
  using namespace std;

  template <class T, class Seq>
  SymMatrix<T> parseSymMatrix(Seq W, std::size_t n) {
    SymMatrix<T> matrix = SymMatrix<T>(n);
    parlay::parallel_for(0, n, [&](size_t i){
        parlay::parallel_for(i+1, n,[&] (size_t j){
	        // matrix.update(i, j, (T)atof(W[i*n + j]));
          matrix.update(i, j, (T)W[i*n + j]);
        });
    });
    return matrix;
  }

  // read a symmatric matrix from file
  template <class T>
  SymMatrix<T> readSymMatrixFromFile(char const *fname, std::size_t n) {
    parlay::sequence<double> W = parlay::sequence<double>(n*n);
    std::ifstream myFile (fname, ios::in | ios::binary);
    myFile.read((char*)W.data(), sizeof(double) * n*n);
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
}  // namespace in_memory
}  // namespace research_graph