#include "dbht.h"

#include <queue>

using namespace std;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::bfs(std::size_t bb, bool *_visited, std::size_t nc, std::size_t* map){
    vtx t1, t2, t3,t4;
    tie(t1,t2,t3,t4) = cliques[bb];

    bool* visited = _visited+bb*nb;

    queue<std::size_t> q;
    q.push(bb);
    while(!q.empty()){
        std::size_t s = q.size();
        while(s>0){
            std::size_t cbb = q.front();
            q.pop();
            visited[cbb] = true;
            
            if(converging[cbb]){
                std::size_t cvgbb = map[cbb];
                
                v2cl[t1*nc + cvgbb] = true;
                v2cl[t2*nc + cvgbb] = true;
                v2cl[t3*nc + cvgbb] = true;
                v2cl[t4*nc + cvgbb] = true;
                b2cl[bb*nc + cvgbb] = true;
            }

            //check parent
            if(cbb!=root){
                std::size_t nextbb = tree[cbb].parent; 
                //look at direction and go to parent if current bb's direction is true
                if(!visited[nextbb] && directions[cbb]){
                    q.push(nextbb);
                }
            }
            //check children
            for(int k=0;k<3;++k){
                if(!tree[cbb].isEmpty(k)){
                    std::size_t nextbb = tree[cbb].children[k]; 
                    //look at direction and go to parent if nextbb's direction is false
                    if(!visited[nextbb] && (!directions[nextbb]))q.push(nextbb);
                }else{
                    break;
                }
            }
            s--;
        }
    }
}


template<class T, class PROF> 
void DBHTTMFG<T, PROF>::nonDiscreteClustering(){
        auto tmp_lst = parlay::delayed_seq<std::size_t>(nb, [&](size_t i){ return i; });
        cvg_bubbles = parlay::filter(tmp_lst, [&](auto i){ 
            return converging[i]; 
        });
        parlay::sort_inplace(cvg_bubbles);
        nc = cvg_bubbles.size();
        parlay::parallel_for(0,nc, [&](std::size_t i){
            map[cvg_bubbles[i]] = i;
        });


        v2cl.resize(n*nc, NO_REACH);
        b2cl.resize(nb*nc, NO_REACH);
        parlay::sequence<bool> visited(nb * nb, false);// helper array for bfs

        parlay::parallel_for(0,nb, [&](std::size_t i){
            bfs(i, visited.begin(), nc, map.begin());
        });
    }
