#include "partmfg.h"
#include "dbht.h"

template<class T, class PROF> 
void ParTMFG<T, PROF>::initGainArray(){
        if(!hasUninsertedV()) return;
        if(use_heap){
            initGainArrayHeap();
            return;
        }
        parlay::parallel_for(0, triangles_ind, [&](face i) {
            auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
            pair<vtx, T> result = getBestGain(in, triangles[i]);
            vtx v = result.first;
            max_clique_gains[i] = make_tuple(v, result.second, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
            face_store[num_faces-1+max_face_num*v] = i;
        });
    }

template<class T, class PROF> 
parlay::sequence<size_t> ParTMFG<T, PROF>::getAllBestVertices(size_t THRESHOLD){
        // get the best triangle for each conflict vertex
        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, THRESHOLD));

        auto tmp_lst = parlay::delayed_seq<size_t>(THRESHOLD, [&](size_t i){ return i; });
        auto insert_list = parlay::filter(tmp_lst, [&](auto i){ 
        return i==THRESHOLD-1 || (get<0>(max_clique_gains[i])!= get<0>(max_clique_gains[i+1])); 
        });
        return insert_list;

    }

template<class T, class PROF> 
parlay::sequence<size_t> ParTMFG<T, PROF>::getBestVertices(size_t THRESHOLD){
#ifdef DEBUG
        assert(THRESHOLD <= triangles_ind);
        checkGains(0);
#endif  
        if(THRESHOLD < triangles_ind){
            parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind), 
            [&](const gainT &i, const gainT &j){
                if(get<1>(i) == get<1>(j)) return i < j;
                return get<1>(i) > get<1>(j); // the largest at the front after sorting
            });
        }
        return getAllBestVertices(THRESHOLD);

    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateGainArray(sequence<size_t> &insert_list){
    if(!hasUninsertedV()) return;
    if(use_heap){
        updateGainArrayHeap(insert_list);
        return;
    }
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
            vtx old_v = v_list[j];
            parlay::parallel_for(0, vtx_to_face_inds[old_v], [&](size_t ii) {
                face i = face_store[ii+max_face_num*old_v];
                pair<vtx, T> result = getBestGain(in, triangles[i]);
                vtx new_v = result.first;
                max_clique_gains[i] = make_tuple(new_v, result.second, i);
                size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
                face_store[num_faces-1+max_face_num*new_v] = i; //TODO: CHANGE DUE TO CONTENTION?
            });
        });
    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateGainArrayHeap(sequence<size_t> &insert_list){
        // auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        auto v_list = sequence<pair<vtx,face>>::uninitialized(insert_list.size());
        parlay::parallel_for(0, insert_list.size(), [&](size_t i) {
            auto entry = max_clique_gains[insert_list[i]];
            v_list[i]= make_pair(get<0>(entry), get<2>(entry));
        });
// optimize this
        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind-2*insert_list.size()), 
            [&](const gainT &i, const gainT &j){
                return get<2>(i) < get<2>(j); // sort by face indices to get the corresponding orders again, sort only the before the insert part
            });
#ifdef DEBUG
        checkGains(2*insert_list.size());
#endif
        size_t triangles_ind_prev = triangles_ind - 2*insert_list.size();
        parlay::parallel_for(0, v_list.size(), [&](size_t j) {
            vtx old_v = v_list[j].first; face tri = v_list[j].second;
            parlay::parallel_for(0, vtx_to_face_inds[old_v], [&](size_t ii) {
                face i = face_store[ii+max_face_num*old_v];
                if(i == tri || i >= triangles_ind_prev ){heapifyFace(i);}
                heapEle result = getMinValidHeapEle(i);
                vtx new_v = result.second;
                max_clique_gains[i] = make_tuple(new_v, result.first, i);
                size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
                face_store[num_faces-1+max_face_num*new_v] = i; 
            });
        });
    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::inertMultiple(sequence<size_t> &insert_list, DBHTTMFG<T, PROF> *clusterer){
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
            cliques[peo_ind + i-3] = cliqueT(t1,t2,t3,v);
            insertToP(v, t1, P_ind + i*3);
            insertToP(v, t2, P_ind + i*3 + 1);
            insertToP(v, t3, P_ind + i*3 + 2);
            vertex_flag[v] = false;
            triangles[tri] = triT(t1,t2,v);
            triangles[triangles_ind + i*2] = triT(t2,t3,v);
            triangles[triangles_ind + i*2 + 1] = triT(t1,t3,v);

            size_t num_faces = vtx_to_face_inds[v];
            face_store[num_faces+max_face_num*v] = triangles_ind + i*2;
            face_store[num_faces+1+max_face_num*v] = triangles_ind + i*2 + 1;
            vtx_to_face_inds[v]+=2;
            if(clusterer != nullptr){
                clusterer->updateDegrees(t1,t2,t3,v);
                clusterer->insertOne(tri, triangles_ind + i*2, triangles_ind + i*2 + 1, peo_ind+i-3, v);
            }
        });
        peo_ind += insert_list.size();
        triangles_ind += 2*insert_list.size();
        P_ind += 3*insert_list.size();
        filterVbyFlag();
    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::insertOne(DBHTTMFG<T, PROF> *clusterer){ // = nullptr
         // get the best vertex among all candidates
        auto *entry_pointer = parlay::max_element(make_slice(max_clique_gains).cut(0,triangles_ind), 
                [&](const auto &i, const auto &j){ 
                    if(get<1>(i) == get<1>(j)) return i > j;
                    // might not be exactly the same as MATLAB, 
                    // because triangle id can be different
                    return get<1>(i) < get<1>(j); 
                });
        auto entry = *entry_pointer;
        vtx v = get<0>(entry);
        face tri = get<2>(entry);
        vtx t1,t2,t3;
        tie(t1,t2,t3) = triangles[tri];

        peo[peo_ind ] = v;
        cliques[peo_ind-3] = cliqueT(t1,t2,t3,v);
        insertToP(v, t1, P_ind);
        insertToP(v, t2, P_ind + 1);
        insertToP(v, t3, P_ind + 2);
        triangles[tri] = triT(t1,t2,v);
        triangles[triangles_ind ] = triT(t2,t3,v);
        triangles[triangles_ind + 1] = triT(t1,t3,v);
        size_t num_faces = vtx_to_face_inds[v];
        face_store[num_faces+max_face_num*v] = triangles_ind;
        face_store[num_faces+1+max_face_num*v] = triangles_ind + 1;

        if(clusterer != nullptr){
            clusterer->updateDegrees(t1,t2,t3,v);
            clusterer->insertOne(tri, triangles_ind, triangles_ind + 1, peo_ind-3, v);
        }

        peo_ind += 1;
        triangles_ind += 2;
        P_ind += 3;
        vtx_to_face_inds[v]+=2;

        removeOneV(v);

        if(vertex_num == 0) return;
if(use_heap){
        vertex_flag[v] = false;
        // get the best vertex and gain for each new triangle
        parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
            face i = face_store[ii+max_face_num*v];
            if(i == tri || i >= triangles_ind -2 ){heapifyFace(i);}
            heapEle result = getMinValidHeapEle(i);
            vtx v = result.second;
            max_clique_gains[i] = make_tuple(v, result.first, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
            face_store[num_faces-1+max_face_num*v] = i;
        });
}else{
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        // get the best vertex and gain for each new triangle
        parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
            face i = face_store[ii+max_face_num*v];
            pair<vtx, T> result = getBestGain(in, triangles[i]);
            vtx v = result.first;
            max_clique_gains[i] = make_tuple(v, result.second, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
            face_store[num_faces-1+max_face_num*v] = i;
        });
} //end else use_heap
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::init(){
        P = sequence<tuple<vtx,vtx,T>>::uninitialized(3*n);
        max_clique_gains = sequence<gainT>::uninitialized(3*n - 6); // ( vertex,gain, face) tuple<vtx, T, face>;
        max_face_num = 3*n - 6;
        cliques = sequence<cliqueT>::uninitialized(n);
        triangles = sequence<triT>::uninitialized(3*n);
        vertex_list = sequence<vtx>::uninitialized(2 *n);
        peo = sequence<vtx>::uninitialized(n);

        vertex_flag = parlay::sequence<bool>(n, true);
        face_store = sequence<size_t>::uninitialized(max_face_num*n);
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

        if(use_heap) initHeap();
        if(use_sorted_list) initSortedList();
    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::initHeap(){
    heap_buffer = sequence<heapEle>::uninitialized(3*n*n); 
    heap_LR = sequence<size_t>::uninitialized(3*n*n);
    heaps=sequence<heapT>::uninitialized(3*n);
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::initSortedList(){
    // heap_buffer = sequence<heapEle>::uninitialized(3*n*n); 
    // heap_LR = sequence<size_t>::uninitialized(3*n*n);
    fv_heaps=sequence<sequence<heapEle>>::uninitialized(n);
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::initGainArrayHeap(){
    parlay::parallel_for(0, triangles_ind, [&](face i) {
        heapifyFace(i);
        heapEle result = getMinValidHeapEle(i);
        vtx v = result.second;
        max_clique_gains[i] = make_tuple(v, result.first, i);
        size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
        face_store[num_faces-1+max_face_num*v] = i;
    });
}

template<class T, class PROF> 
cliqueT ParTMFG<T, PROF>::maxClique(){
        auto Wseq = parlay::delayed_seq<T>(n*n, [&](size_t i) {
            // if(i/n == i%n) return 0.0;
            return getW(i/n, i%n);
        });
        double W_mean = parlay::reduce(Wseq) / n / n;
        // auto W_larger = parlay::delayed_seq<T>(n*n, [&](long i){ return getW(i/n, i%n) > W_mean? getW(i/n, i%n):0; });
        auto W_larger = parlay::delayed_seq<T>(n*n, [&](long i){ return Wseq[i] > W_mean? Wseq[i]:0; });


        sequence<pair<T, vtx>> v_weight = sequence<pair<T, vtx>>(n);
        parlay::parallel_for(0, n, [&](size_t i) {
            int start = i*n; int end = (i+1)*n;
            v_weight[i] = make_pair(parlay::reduce(make_slice(W_larger).cut(start,end)), i);
        });

        parlay::sort_inplace(v_weight); // can improve to 4 max operations

        return cliqueT(get<1>(v_weight[n-1]), get<1>(v_weight[n-2]), get<1>(v_weight[n-3]), get<1>(v_weight[n-4]));
    }