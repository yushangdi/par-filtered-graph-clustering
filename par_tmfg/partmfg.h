#pragma once

#include <tuple>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <algorithm>


#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/delayed_sequence.h"

#include "gettime.h"
#include "atomics.h"
#include "utilities.h"

#include "heap.h"
#include "dbht.h"

using namespace std;
using parlay::sequence;	
using parlay::slice;	

template<class T, class PROF> 
struct ParTMFG{
    using gainT=tuple<vtx, T, face>; 
    using heapEle=pair<T, vtx>;
    using heapT=binary_min_heap<heapEle>;
    // struct entry {
    // using key_t = heapEle;
    // static bool comp(key_t a, key_t b) { 
    //     return a < b;}
    // };
    // using heapT=pam_set<entry>;
    SymM<T>* W;
    sequence<tuple<vtx,vtx,T>> P;
    sequence<gainT> max_clique_gains;

    sequence<cliqueT> cliques;
    sequence<triT> triangles;
    sequence<vtx> peo;
    size_t n;

    PROF *pf;

    ParTMFG(SymM<T> *W_, size_t n_, PROF *_profiler, bool _use_heap=false):W(W_), n(n_), pf(_profiler), use_heap(_use_heap) {
    }

    // initialize the four cliques
    void init();

    //require max_clique_gain updated
    void insertOne(DBHTTMFG<T, PROF> *clusterer=nullptr);
    void initGainArray();


    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    parlay::sequence<size_t> getAllBestVertices(size_t THRESHOLD);

    // get the best vertices of the largest P gains in max_clique_gain vector
    // require: THRESHOLD must be not larger than triangles_ind
    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    parlay::sequence<size_t> getBestVertices(size_t THRESHOLD);

    // update the max_clique_gain array after calling 
    // update the faces in [0, vtx_to_face_inds[v]] of face_store[v*n] for all v in insert_list
    // max_clique_gain does not neccessarily correspond to the order in triangles when start, but does after
    void updateGainArray(sequence<size_t> &insert_list);

    // insert all vertices in max_clique_gain[insert_list].
    // i.e. the indices in insert_list are indices of max_clique_gain
    void inertMultiple(sequence<size_t> &insert_list, DBHTTMFG<T, PROF> *clusterer=nullptr);

    inline T computeGain(vtx i, triT triangle){
        return getW(i,get<0>(triangle)) + getW(i,get<1>(triangle)) + getW(i,get<2>(triangle));
    }

    //get the best vertex for triangle in vertex_list
    pair<vtx, T> getBestGain(parlay::slice<vtx *, vtx *> &vertex_list, triT triangle){
// #ifdef PROFILE
//         timer t;t.start();
// #endif  
        vtx *best_vertex = parlay::max_element(vertex_list, 
            [&](const auto &i, const auto &j){ 
                T ga = computeGain(i,triangle);
                T gb = computeGain(j,triangle); 
                if(ga == gb){
                    return i < j;
                    //want the smallest id vertex if gain equal in matlab, i > j
                }
                return  ga < gb;
            });
        T gain = computeGain(*best_vertex,triangle);
        return make_pair(*best_vertex, gain);
    }

    inline void insertToP(vtx i, vtx j, size_t ind){
        P[ind] = make_tuple(i, j, getW(i, j));
    }

    inline bool hasUninsertedV(){return vertex_num > 0;}    
    inline T getW(vtx i, vtx j){
// #ifdef DEBUG
//         if(i==j){
//             cout << "i==j in getW" << endl;
//         }
// #endif
        //{return W[i*n + j];}
        return W->get(i,j);
    }
    inline size_t getTrianglesNum(){return triangles_ind;}

    T computeCost(){
      // compute the total gain here
        T total_cost = 0;
        for(size_t i=0; i<P_ind; ++ i){
            // print(P[i]);
            total_cost += get<2>(P[i]);
        }

        cout << std::setprecision(20) << 2*total_cost << endl;
        cout << P_ind << endl; 
        return 2*total_cost;   
    }

    void outputPeo(string filename, size_t _n=0){
        if(_n==0)_n = n;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << peo[i]+1 << endl;
        }
        file_obj.close();
    }    

    void outputP(string filename, size_t _n=0){
        if(_n==0)_n = P_ind;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << get<0>(P[i])+1 << " " << get<1>(P[i])+1 << " " << get<2>(P[i]) << endl;
        }
        file_obj << n << " " << n << " " << 0 << endl;
        file_obj.close();
    }   

    void outputCliques(string filename, size_t _n=0){
        if(_n==0)_n = n-3;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << get<0>(cliques[i]) << " " << get<1>(cliques[i]) << " " << get<2>(cliques[i]) << " " << get<3>(cliques[i]) << endl;
        }
        file_obj.close();
    }    
#ifdef DEBUG
        bool validV(vtx t){
            return t < n && t >= 0;
        }

        void checkTriangles(){
            for (size_t i=0; i < triangles_ind; ++i){
                vtx t1,t2,t3;
                tie(t1,t2,t3) = triangles[i];
                if(t1 == t2 || t1==t3 || t2==t3){
                    cout << t1 << " " << t2 << " " << t3 << endl;
                }
            }


        }

        void checkGains(size_t k){
            for (size_t i=0; i < triangles_ind-k; ++i){
                auto entry = max_clique_gains[i];
                face tri = get<2>(entry);
                if(tri != i){
                    cout << tri << " " << i << endl;
                }
            }


        }
#endif

    private:
    parlay::sequence<bool> vertex_flag;
    sequence<vtx> vertex_list;
    sequence<size_t> face_store;
    // sequence<T> W2;

    sequence<size_t> vtx_to_face_inds; //how many faces have vtx as the best vtx
    size_t vertex_num;
    size_t vertex_start;
    size_t triangles_ind;
    size_t peo_ind;
    size_t P_ind = 0;
    size_t max_face_num;
    //// used for heap
    bool use_heap;
    sequence<heapT> heaps;
    sequence<heapEle> heap_buffer;
    sequence<size_t> heap_LR;
    bool use_sorted_list = true; // change to input
    sequence<size_t> sorted_list_pointer; // stores the next vertex to look at


    //allocate space for heap
    void initHeap();
    void initGainArrayHeap();
    void updateGainArrayHeap(sequence<size_t> &insert_list);

    inline heapEle negateGain(heapEle ele){
        return heapEle(-1.0*ele.first, ele.second);
    }

    //prereq: vertex_falg is updated
    inline heapEle getMinValidHeapEle(face i){
// #ifdef PROFILE
//         timer t;t.start();
// #endif  
        pf->incFindMinNum();
        heapEle ele;
        do {
// #ifdef DEBUG
//     if(heaps[i].size == 0){
//         cout << "0 size heap!" <<endl;
//     }
// #endif
            if(use_sorted_list){ 
                size_t ind = sorted_list_pointer[i];
                sorted_list_pointer[i]++;
                ele = heap_buffer[i*n+ind];
            }else{
            size_t ind = heaps[i].argmin();
            ele = heap_buffer[i*n+ind];
            heaps[i].heap_pop();
            // ele = *heaps[i].select(0);
            // heaps[i]=heapT::remove(move(heaps[i]), ele);
            }
        }
        while ( !vertex_flag[ele.second] );
// #ifdef PROFILE
//         // pf->incFindMinTime(t.next());
//         cout << "find min: "<<t.next()<<endl;
// #endif  
        return negateGain(ele);
    }
    
    // init a heap buffer for triangles[i]
    inline void heapifyFace(face i){
// #ifdef PROFILE
//         timer t1;t1.start();
// #endif  
        pf->incHeapifyNum();
        pf->incHeapifySize(vertex_num);
        triT t = triangles[i];
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        parlay::parallel_for(0, in.size(), [&](size_t v_ind) {
            vtx v = in[v_ind];
            T gain = computeGain(v, t);
            heap_buffer[i*n+v_ind] = heapEle(-1.0*gain, v);
        });
        if(use_sorted_list){
        // parlay::sort_inplace(make_slice(heap_buffer).cut(i*n, i*n+vertex_num));
        // parlay::internal::seq_sort_inplace(make_slice(heap_buffer).cut(i*n, i*n+vertex_num), std::less<heapEle>{}, false);
        parlay::internal::quicksort(make_slice(heap_buffer).cut(i*n, i*n+vertex_num), std::less<heapEle>{});
        sorted_list_pointer[i] = 0;
       }else{
        heaps[i] = binary_min_heap<heapEle>(heap_buffer.data()+ (i*n), vertex_num, heap_LR.data()+ (i*n));
        heaps[i].heapify();
        // heaps[i] = heapT();
        // heaps[i] = heapT::multi_insert(move(heaps[i]), make_slice(heap_buffer).cut(i*n, i*n+vertex_num));
       }

// #ifdef PROFILE
//         // pf->incHeapifyTime(t1.next());
//         cout << vertex_num << endl;
//         cout << "heapify: "<<t1.next() << endl;;
// #endif  
    }

    //// used for heap

    //tested
    cliqueT maxClique();

    void removeOneV(vtx v){
        //update vertex list
        //filter and then copy back
        size_t new_vertex_start = 0;
        if(vertex_start == 0){
            new_vertex_start = n;
        }
        auto out2 = make_slice(vertex_list).cut(new_vertex_start, n);
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);

        vertex_num = parlay::filter_into(in, out2, [&](vtx i){ 
            return i != v; 
            });
        vertex_start = new_vertex_start;
    }

    void filterVbyFlag(){
        //update vertex list
        //filter and then copy back
        size_t new_vertex_start = 0;
        if(vertex_start == 0){
            new_vertex_start = n;
        }
        auto out2 = make_slice(vertex_list).cut(new_vertex_start, n);
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);

        vertex_num = parlay::filter_into(in, out2, [&](vtx i){ 
            return vertex_flag[i]; 
            });
        vertex_start = new_vertex_start;
    }
};

