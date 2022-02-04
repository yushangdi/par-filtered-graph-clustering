#include "dbht.h"

// #include <pair>
using namespace std;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::assignToConvergingBubble(){
        // if(nc==1){ // if there is one converging bubble, all vertices belong to a single cluster
        //     return;
        // }

        assigned.resize(n, false);
        //store the converging bubble membership of the vtx
        AT* vtxInConv = (AT*)malloc(n * sizeof(AT));
         parlay::parallel_for(0, n, [&](std::size_t i) {
            vtxInConv[i] = AT(NO_ATTACH_MIN, SIZE_T_MAX);
         });

        parlay::parallel_for(0, n-3, [&](std::size_t i) {
            // compute the distance to this convergin bubble
            // if smaller, update
            if(converging[i]){
                vtx u = get<0>(cliques[i]);
                T chi = computeChi(u, i);
                pbbs::write_max(&vtxInConv[u], AT(chi,i), cAT());

                u = get<1>(cliques[i]);chi = computeChi(u, i);
                pbbs::write_max(&vtxInConv[u], AT(chi,i), cAT());

                u = get<2>(cliques[i]);chi = computeChi(u, i);
                pbbs::write_max(&vtxInConv[u], AT(chi,i), cAT());

                u = get<3>(cliques[i]);chi = computeChi(u, i);
                pbbs::write_max(&vtxInConv[u], AT(chi,i), cAT());
            }
        });

        // sort the vtx By bubble membership
        auto vtxInConvSort = parlay::sequence<pair<std::size_t, vtx>>(n); // bbid, v
        auto offset = parlay::sequence<std::size_t>(nc+1, n);
         parlay::parallel_for(0, n, [&](vtx i) {
            std::size_t id = vtxInConv[i].second;
            vtxInConvSort[i] = make_pair(id, i);
            bbMember[i] = id; // should store the bubble id, not the cluster id
            if(id != SIZE_T_MAX) assigned[i] = true;
         });
        parlay::sort_inplace(vtxInConvSort);
        offset[map[vtxInConvSort[0].first]] = 0;
        parlay::parallel_for(0, n-1, [&](vtx i) {
            if(vtxInConvSort[i].first!=vtxInConvSort[i+1].first){
                std::size_t id = vtxInConvSort[i+1].first;
                if(id==SIZE_T_MAX){
                    id = nc;
                }else{
                    id = map[id];
                }
                offset[id] = i+1;
            }
        });
        parlay::parallel_for(0, nc-1, [&](std::size_t i) {
            if(offset[i]==n){
                offset[i]=offset[i+1];
            }
        });
        
        //Compute the distance between a vertex and the converging bubbles.
        //Look for the closest converging bubble
        //Assign discrete cluster membership according to the distances to the converging bubbles
        // in the matlab code, if there is an empty converging bubble, the vertices that can only reach
        // this converging bubble will not be ablt to be assigned to any cluster, 
        // and there will be an error in the program
         

        parlay::parallel_for(offset[nc], n, [&](vtx i) { //have not been assigned to any converging bubble
            std::size_t v = (std::size_t)vtxInConvSort[i].second;
            vtxInConv[v] = AT(NO_ATTACH_MAX, SIZE_T_MAX);  //RESET the -1 ones to max for write min
            parlay::parallel_for(0, nc, [&](std::size_t j) { // compute the mean distance from remainning vertices to all converging bubbles
                if(v2cl[v*nc + j]!= NO_REACH){
                    T avg = 0;
                    if((offset[j+1] - offset[j]) > 0){ // not an empty cluster
                        for(std::size_t k = offset[j];k< offset[j+1];++k){
                            avg += SP.get(v, vtxInConvSort[k].second);
                        }
                        avg /= (offset[j+1] - offset[j]);
                    }else{
                        avg = NO_ATTACH_MAX;
                    }
                    pbbs::write_min(&vtxInConv[v], AT(avg,cvg_bubbles[j]), cAT());
                }
            });
            // if(vtxInConv[v].first == NO_ATTACH_MAX){
            //     cout << "v not assigned " << v << " " << vtxInConv[v].first << " " << vtxInConv[v].second<< endl;
            // }
        });

         parlay::parallel_for(0, n, [&](vtx i) {
            // flatClustering[i] = vtxInConv[i].load().second;
            flatClustering[i] = vtxInConv[i].second;
         });

        free(vtxInConv);

    }

// the paper is different from their implmenetation
// the paper says divided by the size of the bubble
// but the MATLAB implementation divided by the sum of the edge weights in the bubble
// also, all vertices' bubble membership need to be computed, this
// is different from the paper, where vertices assigned to V^0 are not re-assigned
// also, does not only assign to bubbles that can reach the converging bubble (in the subtree)
// moreover, the max might be resolved to be a different bubble if two have the same attachment
template<class T, class PROF> 
void DBHTTMFG<T, PROF>::assignToBubble(){

    AT* vtxInB = (AT*)malloc(n * sizeof(AT));
    parlay::parallel_for(0, n, [&](std::size_t i) {
        vtxInB[i] = AT(NO_ATTACH_MIN, SIZE_T_MAX);
    });

    parlay::parallel_for(0,nb,[&](std::size_t bb){
        T total_chi = 0; // should be divided by 2 in the end, omit because we only care about max
        vtx v0,v1, v2, v3; tie(v0,v1,v2,v3) = cliques[bb];

        T chi_raw0 = computeChi(v0, bb); total_chi += chi_raw0;
        T chi_raw1 = computeChi(v1, bb); total_chi += chi_raw1;
        T chi_raw2 = computeChi(v2, bb); total_chi += chi_raw2;
        T chi_raw3 = computeChi(v3, bb); total_chi += chi_raw3;

        if(true){  pbbs::write_max(&vtxInB[v0], AT(chi_raw0/total_chi,bb), cAT());}
        if(true){  pbbs::write_max(&vtxInB[v1], AT(chi_raw1/total_chi,bb), cAT());}
        if(true){  pbbs::write_max(&vtxInB[v2], AT(chi_raw2/total_chi,bb), cAT());}
        if(true){  pbbs::write_max(&vtxInB[v3], AT(chi_raw3/total_chi,bb), cAT());}
    });

    parlay::parallel_for(0, n, [&](vtx i) {
        bbMember[i] = vtxInB[i].second;
    });

    free(vtxInB);
}

// version where vertices in converging bubbles keep their assignments
// and only assign to bubbles that can reach the converging bubble (in the subtree)
// template<class T, class PROF> 
// void DBHTTMFG<T, PROF>::assignToBubble(){

//     auto vtxInB = parlay::sequence<atomic<AT>>(n);
//     parlay::parallel_for(0, n, [&](std::size_t i) {
//         if(!assigned[i])vtxInB[i].store(AT(NO_ATTACH_MIN, SIZE_T_MAX));
//     });

//     parlay::parallel_for(0,nb,[&](std::size_t bb){
//         T total_chi = 0; // should be divided by 2 in the end, omit because we only care about max
//         vtx v0,v1, v2, v3; tie(v0,v1,v2,v3) = cliques[bb];
//         T chi_raw0 = computeChi(v0, bb); total_chi += chi_raw0;
//         T chi_raw1 = computeChi(v1, bb); total_chi += chi_raw1;
//         T chi_raw2 = computeChi(v2, bb); total_chi += chi_raw2;
//         T chi_raw3 = computeChi(v3, bb); total_chi += chi_raw3;

//         if((!assigned[v0]) && b2cl[bb*nc + map[flatClustering[v0]]]){  pbbs::write_max(&vtxInB[v0], AT(chi_raw0/total_chi,bb), cAT());}
//         if((!assigned[v1]) && b2cl[bb*nc + map[flatClustering[v1]]]){  pbbs::write_max(&vtxInB[v1], AT(chi_raw1/total_chi,bb), cAT());}
//         if((!assigned[v2]) && b2cl[bb*nc + map[flatClustering[v2]]]){  pbbs::write_max(&vtxInB[v2], AT(chi_raw2/total_chi,bb), cAT());}
//         if((!assigned[v3]) && b2cl[bb*nc + map[flatClustering[v3]]]){  pbbs::write_max(&vtxInB[v3], AT(chi_raw3/total_chi,bb), cAT());}
//     });

//     parlay::parallel_for(0, n, [&](vtx i) {
//        if(!assigned[i])bbMember[i] = vtxInB[i].load().second;
//     });
// }