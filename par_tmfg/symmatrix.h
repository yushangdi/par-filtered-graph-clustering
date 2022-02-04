#pragma once

using namespace std;

// symmetric matrix, diagonal is all diag_val = 0 by default
template<class T> 
struct SymM{

    std::size_t n = 0;
    T* distmat = nullptr;
    std::size_t distmat_size = 0;
    T diag_val = 0;

    SymM(){}
    SymM(std::size_t _n, T _diag_val = 0){
        n = _n;
        diag_val = _diag_val;
        distmat_size = n*(n-1)/2;
        distmat = (T *)malloc(distmat_size * sizeof(T));
    }

    void init(std::size_t _n){
        n = _n;
        distmat_size = n*(n-1)/2;
        distmat = (T *)malloc(distmat_size * sizeof(T));
    }

    inline long getInd(std::size_t i, std::size_t j){
        if(i == j) return distmat_size;
        long r_ = static_cast<long>(i);
        long c_ = static_cast<long>(j);
        return (((2*n-r_-3) * r_) >> 1 )+ (c_)-1;
    }
    
    inline void setDiag(T v){
        diag_val = v;
    }

    inline T get(std::size_t r_, std::size_t c_){
        if(r_ == c_) return diag_val;
        if(r_ > c_) swap(r_,c_);
        return( distmat[getInd(r_, c_)] );
    }

    inline void update(std::size_t r_, std::size_t c_, T dist){
        if(r_ == c_) return;
        if(r_ > c_) swap(r_,c_);
        distmat[getInd(r_, c_)] = dist;
    }

    ~SymM(){
        free(distmat);
    }
    void printMatrix(){
    // for (std::size_t i=0; i<distmat_size; i++){
    //     cout << distmat[i] << endl;
    // }
    for (std::size_t i=0; i<n; i++){
        for (std::size_t j=0; j<n; j++) {
            cout << get(i,j) << " ";
        }
        cout << endl;
    } 
    cout << "===" << endl;
    for (std::size_t i=0; i<n; i++){
        for (std::size_t j=i+1; j<n; j++) {
            cout << i << " " << j << " " <<getInd(i,j) << " " << get(i,j) << endl;
        }
    } 
    }
};
