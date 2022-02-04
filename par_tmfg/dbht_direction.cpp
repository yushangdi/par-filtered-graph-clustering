#include "dbht.h"

using namespace std;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::computeDirection(){
        // computeDirectionHelper(root);
        parlay::parallel_for(0, 3, [&](int k) {
            if(!tree[root].isEmpty(k)){
                std::size_t cid = tree[root].children[k];
                computeDirectionHelper(cid, root);
            }
        }, 1);
    }


template<class T, class PROF> 
typename DBHTTMFG<T, PROF>::partialDeg DBHTTMFG<T, PROF>::computeDirectionHelper(std::size_t id, std::size_t pid){
        partialDeg result;

        vtx t1,t2,t3,v; tie(t1,t2,t3,v) = compareCliques(id, pid);
        result.inds[0] = t1; result.inds[1] = t2; result.inds[2] = t3;
        result.vals[0] += getW(t1, v);
        result.vals[1] += getW(t2, v);
        result.vals[2] += getW(t3, v);

        parlay::parallel_for(0, 3, [&](int k) {
            if(!tree[id].isEmpty(k)){
                std::size_t cid = tree[id].children[k];
                partialDeg cresult = computeDirectionHelper(cid, id);
                result.inc(cresult);
            }
        }, 1);
    
        if(id != root){
            T tri_val = getW(t1, t2) + getW(t1, t3) + getW(t2, t3);
            T in_val = result.vals[0] + result.vals[1] + result.vals[2]; 
            T out_val = degrees[t1] +  degrees[t2] +  degrees[t3] - in_val - 2*tri_val;
            if(in_val > out_val){ // goinig from parent to this bubble
                directions[id] = false;
                std::size_t pid = tree[id].parent;
                if(converging[pid]) converging[pid] = false;
            }else{
                directions[id] = true;
                if(converging[id]) converging[id] = false;
            }
        }
        
        return result;//make_tuple(vals[0], vals[1], vals[2]);
    }

template<class T, class PROF> 
cliqueT DBHTTMFG<T, PROF>::compareCliques(std::size_t id1, std::size_t id2){
        vtx inds3[4];
        int ct = 0;
#ifdef DEBUG
        bool flag = false;
#endif
        if(cChelper<0>(id1, id2)){
            inds3[ct]=get<0>(cliques[id1]);ct++;
        }else{
            inds3[3]=get<0>(cliques[id1]);
        }

        if(cChelper<1>(id1, id2)){
            inds3[ct]=get<1>(cliques[id1]);ct++;
        }else{
#ifdef DEBUG
        assert(!flag);
        flag = true;
#endif
            inds3[3]=get<1>(cliques[id1]);
        }

        if(cChelper<2>(id1, id2)){
            inds3[ct]=get<2>(cliques[id1]);ct++;
        }else{
#ifdef DEBUG
        assert(!flag);
        flag = true;
#endif
            inds3[3]=get<2>(cliques[id1]);
        }

        if(cChelper<3>(id1, id2)){
            inds3[ct]=get<3>(cliques[id1]);ct++;
        }else{
#ifdef DEBUG
        assert(!flag);
        flag = true;
#endif
            inds3[3]=get<3>(cliques[id1]);
        }
#ifdef DEBUG
        assert(ct==3);
#endif
        return cliqueT(inds3[0],inds3[1],inds3[2],inds3[3]);

    }