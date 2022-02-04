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

using namespace std;
using parlay::sequence;	
using parlay::slice;	

using vtx=int;
using face=int;
using cliqueT=tuple<vtx,vtx,vtx,vtx>;
using triT=tuple<vtx,vtx,vtx>;


template<class T> 
struct ParTMFG{
    using gainT=tuple<vtx, T, face>; 
    sequence<T> W;
    sequence<tuple<vtx,vtx,T>> P;
    sequence<gainT> max_clique_gains;

    sequence<cliqueT> cliques;
    sequence<triT> triangles;
    sequence<vtx> peo;
    size_t n;

    ParTMFG(std::string filename){
    size_t m;

    // assnme n^2 memory possible because W is n^2 memory
    tie(m,n) = readMatrix<T>(filename, W);

    assert(m==n);
    }

    // initialize the four cliques
    void init(){
        P = sequence<tuple<vtx,vtx,T>>::uninitialized(3*n);
        max_clique_gains = sequence<gainT>::uninitialized(3*n - 6); // ( vertex,gain, face) tuple<vtx, T, face>;
        cliques = sequence<cliqueT>::uninitialized(n);
        triangles = sequence<triT>::uninitialized(3*n);
        vertex_list = sequence<vtx>::uninitialized(2 *n);
        peo = sequence<vtx>::uninitialized(n);

        vertex_flag = parlay::sequence<bool>(n, true);
        face_store = sequence<size_t>::uninitialized(n*n);
        vtx_to_face_inds = sequence<size_t>(n, 0);

        cliques[0] = maxClique();
        vtx t1, t2, t3,t4;
        tie(t1,t2,t3,t4) = cliques[0];
        peo[0] = t1; peo[1] = t2; peo[2] = t3; peo[3] = t4;
        for(int i=0; i<4; ++i){
            vertex_flag[peo[i]] = false;
            for(int j=i+1; j<4; ++j){
                insertToP(peo[i], peo[j], P_ind);
                P_ind++;
            }
        }
        triangles[0] = triT(t1, t2, t3);
        triangles[1] = triT(t1, t2, t4);
        triangles[2] = triT(t1, t3, t4);
        triangles[3] = triT(t2, t3, t4);

        auto all_vertices = parlay::delayed_seq<vtx>(n, [&](size_t i){ return i; });
        vertex_num = parlay::filter_into(all_vertices, vertex_list, [&](vtx i){ 
            return i != t1 && i != t2 && i != t3 && i != t4; 
            });
        vertex_start = 0;
        triangles_ind = 4;
        peo_ind = 4;
    }

    void initGainArray(){
        
        parlay::parallel_for(0, triangles_ind, [&](face i) {
            auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
            pair<vtx, T> result = getBestGain(in, triangles[i]);
            vtx v = result.first;
            max_clique_gains[i] = make_tuple(v, result.second, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
            face_store[num_faces-1+n*v] = i;
        });
    }

    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    auto getAllBestVertices(size_t THRESHOLD){
        // get the best triangle for each conflict vertex
        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, THRESHOLD));

        auto tmp_lst = parlay::delayed_seq<size_t>(THRESHOLD, [&](size_t i){ return i; });
        auto insert_list = parlay::filter(tmp_lst, [&](auto i){ 
        return i==THRESHOLD-1 || (get<0>(max_clique_gains[i])!= get<0>(max_clique_gains[i+1])); 
        });
        return insert_list;

    }

    // get the best vertices of the largest P gains in max_clique_gain vector
    // require: THRESHOLD must be not larger than triangles_ind
    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    auto getBestVertices(size_t THRESHOLD){
#ifdef DEBUG
        assert(THRESHOLD <= triangles_ind);
#endif  
        if(THRESHOLD < triangles_ind){
            parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind), 
            [&](const gainT &i, const gainT &j){
                return get<1>(i) > get<1>(j); // the largest at the front after sorting
            });
        }

        return getAllBestVertices(THRESHOLD);

    }

    // update the max_clique_gain array after calling 
    // update the faces in [0, vtx_to_face_inds[v]] of face_store[v*n] for all v in insert_list
    // max_clique_gain does not neccessarily correspond to the order in triangles when start, but does after
    void updateGainArray(sequence<size_t> &insert_list){
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        auto v_list = sequence<vtx>::uninitialized(insert_list.size());
        parlay::parallel_for(0, insert_list.size(), [&](size_t i) {
            auto entry = max_clique_gains[insert_list[i]];
            vtx v = get<0>(entry);
            v_list[i]= v;
        });

        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind-2*insert_list.size()), 
            [&](const gainT &i, const gainT &j){
                return get<2>(i) < get<2>(j); // sort by face indices to get the corresponding orders again, sort only the before the insert part
            });
#ifdef DEBUG
        checkGains(2*insert_list.size());
#endif
        parlay::parallel_for(0, v_list.size(), [&](size_t j) {
            vtx v = v_list[j];
            parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
                face i = face_store[ii+n*v];
                pair<vtx, T> result = getBestGain(in, triangles[i]);
                vtx v = result.first;
                max_clique_gains[i] = make_tuple(v, result.second, i);
                size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
                face_store[num_faces-1+n*v] = i; //TODO: CHANGE DUE TO CONTENTION?
            });
        });
    }
    // insert all vertices in max_clique_gain[insert_list].
    // i.e. the indices in insert_list are indices of max_clique_gain
    void inertMultiple(sequence<size_t> &insert_list){
        // for each vertex, triangle pair in the insert_list, insert the vertex 
        // and update triangle list, vertex_flag,  and P
        parlay::parallel_for(0, insert_list.size(), [&](size_t i) {
            auto entry = max_clique_gains[insert_list[i]];
            vtx v = get<0>(entry);
            face tri = get<2>(entry);
            vtx t1,t2,t3;
            tie(t1,t2,t3) = triangles[tri];
#ifdef DEBUG
            if(!(validV(v) && validV(t1) && validV(t2) && validV(t3))){
                cout << v << " " << t1 << " " << t2 << " " << t3 << endl;
            }
#endif

            peo[peo_ind + i] = v;
            insertToP(v, t1, P_ind + i*3);
            insertToP(v, t2, P_ind + i*3 + 1);
            insertToP(v, t3, P_ind + i*3 + 2);
            vertex_flag[v] = false;
            triangles[tri] = triT(t1,t2,v);
            triangles[triangles_ind + i*2] = triT(t2,t3,v);
            triangles[triangles_ind + i*2 + 1] = triT(t1,t3,v);

            size_t num_faces = vtx_to_face_inds[v];
            face_store[num_faces+n*v] = triangles_ind + i*2;
            face_store[num_faces+1+n*v] = triangles_ind + i*2 + 1;
            vtx_to_face_inds[v]+=2;
        });
        peo_ind += insert_list.size();
        triangles_ind += 2*insert_list.size();
        P_ind += 3*insert_list.size();
        filterVbyFlag();
    }


    T getW(vtx i, vtx j){
        return W[i*n + j];
    }

    //tested
    cliqueT maxClique(){
        double W_mean = parlay::reduce(W) / n / n;
        auto W_larger = parlay::delayed_seq<T>(n*n, [&](long i){ return W[i] > W_mean? W[i]:0; });

        sequence<pair<T, vtx>> v_weight = sequence<pair<T, vtx>>(n);
        parlay::parallel_for(0, n, [&](size_t i) {
            int start = i*n; int end = (i+1)*n;
            v_weight[i] = make_pair(parlay::reduce(make_slice(W_larger).cut(start,end)), i);
        });

        parlay::sort_inplace(v_weight); // can improve to 4 max operations

        return cliqueT(get<1>(v_weight[n-1]), get<1>(v_weight[n-2]), get<1>(v_weight[n-3]), get<1>(v_weight[n-4]));
    }

    T computeGain(vtx i, triT triangle){
        return getW(i,get<0>(triangle)) + getW(i,get<1>(triangle)) + getW(i,get<2>(triangle));
    }

    //get the best vertex for triangle in vertex_list
    pair<vtx, T> getBestGain(parlay::slice<vtx *, vtx *> &vertex_list, triT triangle){
        vtx *best_vertex = parlay::max_element(vertex_list, 
            [&](const auto &i, const auto &j){ 
                return computeGain(i,triangle) < computeGain(j,triangle); });
        T gain = computeGain(*best_vertex,triangle);
        return make_pair(*best_vertex, gain);
    }

    void insertToP(vtx i, vtx j, size_t ind){
        P[ind] = make_tuple(i, j, getW(i, j));
    }

    bool hasUninsertedV(){
        return vertex_num > 0;
    }

    T computeCost(){
      // compute the total gain here
        T total_cost = 0;
        for(int i=0; i<P_ind; ++ i){
            // print(P[i]);
            total_cost += get<2>(P[i]);
        }

        cout << 2*total_cost << endl;
        cout << P_ind << endl; 
        return 2*total_cost;   
    }

    void outputPeo(string filename){
        ofstream file_obj;
        file_obj.open(filename); //"+ to_string(round) + "
        for(int i=0;i<n;i++){
            file_obj << peo[i]+1 << endl;
        }
    }
        //require max_clique_gain updated
    void insertOne(){
        // get the best vertex among all candidates
        auto *entry_pointer = parlay::max_element(make_slice(max_clique_gains).cut(0,triangles_ind), 
                [&](const auto &i, const auto &j){ return get<1>(i) < get<1>(j); });
        auto entry = *entry_pointer;
        vtx v = get<0>(entry);
        face tri = get<2>(entry);
        vtx t1,t2,t3;
        tie(t1,t2,t3) = triangles[tri];

        peo[peo_ind ] = v;
        insertToP(v, t1, P_ind);
        insertToP(v, t2, P_ind + 1);
        insertToP(v, t3, P_ind + 2);
        triangles[tri] = triT(t1,t2,v);
        triangles[triangles_ind ] = triT(t2,t3,v);
        triangles[triangles_ind + 1] = triT(t1,t3,v);
        size_t num_faces = vtx_to_face_inds[v];
        face_store[num_faces+n*v] = triangles_ind;
        face_store[num_faces+1+n*v] = triangles_ind + 1;
            
        peo_ind += 1;
        triangles_ind += 2;
        P_ind += 3;
        vtx_to_face_inds[v]+=2;

        removeOneV(v);

        if(vertex_num == 0) return;
        // get the best vertex and gain for each new triangle
        parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
            auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
            face i = face_store[ii+n*v];
            pair<vtx, T> result = getBestGain(in, triangles[i]);
            vtx v = result.first;
            max_clique_gains[i] = make_tuple(v, result.second, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
            face_store[num_faces-1+n*v] = i;
        });

    }

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

    size_t getTrianglesNum(){
        return triangles_ind;
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
                    cout << tri << "" << i << endl;
                }
            }


        }
#endif

    private:
    parlay::sequence<bool> vertex_flag;
    sequence<vtx> vertex_list;
    sequence<size_t> face_store;

    sequence<size_t> vtx_to_face_inds; //how many faces have vtx as the best vtx
    size_t vertex_num;
    size_t vertex_start;
    size_t triangles_ind;
    size_t peo_ind;
    size_t P_ind = 0;
};

