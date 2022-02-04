#include "dbht.h"

using namespace std;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::updateDegrees(vtx t1, vtx t2, vtx t3, vtx v){
        T val = getW(t1, v); T valtotal = val;  pbbs::write_add(&degrees[t1], val); // degrees[t1] += getW(t1, v);
        val = getW(t2, v); valtotal += val;     pbbs::write_add(&degrees[t2], val); // degrees[t2] += getW(t2, v);
        val = getW(t3, v); valtotal += val;     pbbs::write_add(&degrees[t3], val); // degrees[t3] += getW(t3, v);
        degrees[v] += valtotal; // pbbs::write_add(&degrees[v], valtotal);
    }

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::insertOne(std::size_t t1, std::size_t t2, std::size_t t3, std::size_t newc, vtx v){
        std::size_t clique_ind = t2c[t1]; // index of the clique being created, v is inserted in triangles[t1]
        t2c[t1] = newc; t2c[t2] = newc; t2c[t3] = newc; // update the clique that the new triangles correspond to
        if(t1 ==0){ // if inserting into the outer face
            root = newc; // the new clique has an outer face
            tree[clique_ind].parent = newc;
            tree[newc].addChild(clique_ind);
        }else{
            tree[newc].parent = clique_ind;
            tree[clique_ind].addChildPar(newc);
        }
    }

