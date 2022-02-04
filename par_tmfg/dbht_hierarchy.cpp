#include "dbht.h"
#include "unionfind.h"
#include "linkage.h"
#include "chain.h"


template<class TF>
inline void chain_find_nn(std::size_t chainNum, TF *finder, TreeChainInfo *info, 
                          parlay::sequence<std::size_t> &neighbor, parlay::sequence<std::size_t> &offsets,
                                parlay::sequence<std::size_t> &sizes){
    // update edges array
  parlay::parallel_for(0, chainNum, [&](std::size_t i){
    std::size_t cid = info->terminal_nodes[i];
    finder->edges[cid] = EDGE(cid,NO_NEIGH,numeric_limits<double>::max());
#ifdef DEBUG
    if(sizes[offsets[cid]]==0){
      cout << "zero size neighbor list" << endl;
    }
#endif
    if(sizes[offsets[cid]] > 1){
        auto c_neigh = make_slice(neighbor).cut(offsets[cid], offsets[cid]+sizes[offsets[cid]]);
        finder->getNN(cid, c_neigh);
    }
  });
  // update chain
  parlay::parallel_for(0, chainNum, [&](std::size_t i){
      std::size_t cid = info->terminal_nodes[i];
      info->updateChain(cid, finder->edges[cid].second, finder->edges[cid].getW());
  });
}


template<class TF>
inline void link_terminal_nodes(TF *finder, TreeChainInfo *info, std::size_t round, parlay::sequence<std::size_t> &flags,
                                STAGE stage){
  std::size_t chainNum = info->chainNum;
  auto uf = finder->uf;
  EDGE *edges = finder->edges.data();

  parlay::parallel_for(0, chainNum, [&](std::size_t i){
    std::size_t cid = info->terminal_nodes[i];
    if(edges[cid].second != NO_NEIGH){ // do not link if there no neighbor in this round
        std::size_t nn = edges[cid].second; //uf->find()
        std::size_t nn_of_nn =info->getNN(nn);
        // if( nn_of_nn> 0) nn_of_nn = uf->find(nn_of_nn); do not need because all clusters on chain did not merge in the last round

        // avoid merging twice
        // update chain info and merge trees
        if ((!(info->isTerminal(nn) && cid <= nn)) && (nn_of_nn == cid)){// if they are RNN //nn != cid && 
        std::size_t newc = uf->link(cid, nn, edges[cid].getW());
        info->invalidate(nn, NO_NEIGH);
        info->invalidate(cid, NO_NEIGH);
        info->invalidateRev(newc);
        std::size_t cluster_num = finder->getNode(cid)->cluster_num + finder->getNode(nn)->cluster_num;
        std::size_t newidx = finder->merge(cid, nn, newc, round, edges[cid].getW());
        finder->nodes[newidx].stage = stage;
        if(stage == STAGE::INTERC){
          finder->nodes[newidx].cluster_num = cluster_num;
          finder->nodes[newidx].stage = stage;
        }
        flags[i] = newc;
        }else{
        flags[i] = NO_NEIGH;
        }
    }else{
       flags[i] = NO_NEIGH; 
    }
  });
  
  auto merged = parlay::filter(make_slice(flags).cut(0, chainNum), [&](std::size_t i){return i!=NO_NEIGH;});

  // insert to new hashtables and delete old hashtables
  parlay::parallel_for(0, merged.size(), [&](std::size_t i){
    std::size_t newc = merged[i];
    finder->updateDist(newc, round);
  });

}


template <PARLAY_RANGE_TYPE R, typename BinaryPredicate, PARLAY_RANGE_TYPE R_out>
auto unique_into(const R& r, BinaryPredicate eq, R_out&& out) {
  auto b = parlay::tabulate(
    parlay::size(r), [&eq, it = std::begin(r)](size_t i)
      { return (i == 0) || !eq(it[i], it[i - 1]); });
  return parlay::pack_into(r, b, out);
}

template<PARLAY_RANGE_TYPE R, PARLAY_RANGE_TYPE R_out>
auto unique_into(const R& r, R_out&& out) {
  return unique_into(r, std::equal_to<parlay::range_value_type_t<R>>(), out);
}

template<class TF>
inline void chain_linkage_matrix(TF *finder, TreeChainInfo* info, 
                                parlay::sequence<std::size_t> &neighbor,
                                parlay::sequence<std::size_t> &offsets,
                                parlay::sequence<std::size_t> &sizes,
                                parlay::sequence<std::size_t> &flags, 
                                STAGE stage, std::size_t &round){
  
  std::size_t n = finder->n;
  std::size_t chainNum = info->chainNum;
  
  auto neighborCopy = parlay::sequence<std::size_t>(n);
  auto chunks = parlay::pack_index(parlay::delayed_seq<bool>(n,[&](std::size_t i){
      if(i==0) return true;
      return (offsets[neighbor[i]]!=offsets[neighbor[i-1]]);
  }));
  std::size_t numDone = chunks.size();

  //pack to be contiguous, do no need to intra-bubble
  if(chainNum < n){
    parlay::parallel_for(0, chunks.size(), [&](std::size_t i){
      std::size_t s = chunks[i];
      if(sizes[s]>1){
        size_t e = s+sizes[s];
        auto out = make_slice(neighborCopy).cut(s, e);
        std::size_t new_size = parlay::filter_into(make_slice(neighbor).cut(s, e), out, 
                    [&](std::size_t i){return finder->isActive(i);});
        e = s+new_size;
        auto out2 = make_slice(neighbor).cut(s, e);
        new_size = unique_into(make_slice(neighborCopy).cut(s, e), out2);
        sizes[s] = new_size;
      }
    });
  }

  while(finder->C > numDone){
    round ++;
    if(round >= n){
        cout << "too many rounds" << endl;
        cout << finder->C << endl;
        exit(1);
    }
    chain_find_nn(chainNum, finder, info, neighbor, offsets, sizes);
    link_terminal_nodes(finder, info, round, flags, stage);
    // get ready for next round
    finder->updateActiveClusters(round);
    info->next(finder);
    chainNum = info->chainNum;

    //pack neighbor and get new offsets
    parlay::copy(neighbor, neighborCopy);
    parlay::parallel_for(0, chunks.size(), [&](std::size_t i){
      std::size_t s = chunks[i];
      if(sizes[s]>1){
        size_t e = s+sizes[s];
        auto out = make_slice(neighbor).cut(s, e);
        std::size_t new_size = parlay::filter_into(make_slice(neighborCopy).cut(s, e), out, 
                    [&](std::size_t i){return finder->isActive(i);});
        sizes[s] = new_size;
      }
    });
  }
}

struct nodeComparatorDBHT{
  double eps;
  std::size_t *flatClustering;
  std::size_t *bbMember;
  nodeComparatorDBHT(double _eps, std::size_t *f, std::size_t *b):eps(_eps), flatClustering(f), bbMember(b) {
  }

  bool operator() (const Node i, const Node j) const {
    pair<std::size_t, STAGE> a = make_pair(flatClustering[i.cId], i.stage);
    pair<std::size_t, STAGE> b = make_pair(flatClustering[j.cId], j.stage);
    if(a!=b){return a < b;}
    if(i.stage != STAGE::INTRABB){ //INTERBB 
      if(abs(i.getHeight() - j.getHeight()) <= eps) return i.getRound() < j.getRound();
      return i.getHeight() < j.getHeight();
    }
    //stage == STAGE::INTRABB 
      if(bbMember[i.cId] != bbMember[j.cId]) return bbMember[i.cId] < bbMember[j.cId];
      if(abs(i.getHeight() - j.getHeight()) <= eps) return i.getRound() < j.getRound();
      return i.getHeight() < j.getHeight();
  }
};

//interClusterID: the first INTERC node
dendroLine* formatDendrogram(parlay::sequence<Node> &nodes, std::size_t n, double eps, std::size_t interClusterID, std::size_t *flatClustering, std::size_t *b){

    // sort by cluster id, within each cluster, sort by bubble id and then height. Intercluster is after intracluster within each cluster
    auto sorted_nodes_front = parlay::sort(make_slice(nodes).cut(n,interClusterID), nodeComparatorDBHT(eps, flatClustering, b));
    auto sorted_nodes_back = parlay::sort(make_slice(nodes).cut(interClusterID,2*n-1), nodeComparator(eps));


    // compute height
    auto chunks = parlay::pack_index(parlay::delayed_seq<bool>(sorted_nodes_front.size(),[&](std::size_t i){
      if(i==0) return true;
      return (flatClustering[sorted_nodes_front[i].cId] != flatClustering[sorted_nodes_front[i-1].cId]);
    }));
    parlay::parallel_for(0, chunks.size(), [&](std::size_t i){
      std::size_t start = chunks[i];
      std::size_t end = i==chunks.size()-1? sorted_nodes_front.size():chunks[i+1];
      std::size_t s = end-start;
      parlay::parallel_for(0, s, [&](std::size_t j){
        sorted_nodes_front[start+j].height = 1.0/(s-j);
      });
    });

    auto sorted_nodes = parlay::delayed_seq<Node>(n-1, [&](std::size_t i){
      if(i < interClusterID-n) return sorted_nodes_front[i];
      return sorted_nodes_back[i-(interClusterID-n)];
    });
    auto map = parlay::sequence<std::size_t>(n);
    parlay::parallel_for(0, n-1, [&](std::size_t i){
        map[sorted_nodes[i].getIdx() - n] = i+n;
    });
    dendroLine* dendrogram = (dendroLine *)malloc(sizeof(dendroLine)*(n-1));//newA(dendroLine, n-1);
    parlay::parallel_for(0, n-1, [&](std::size_t i){
        std::size_t left = sorted_nodes[i].left->getIdx();
        std::size_t right = sorted_nodes[i].right->getIdx();

        left = left < n ? left : map[left-n];
        right = right < n ? right : map[right-n];

        if(left > right) swap(left, right);

        double height = sorted_nodes[i].getHeight();
        if(sorted_nodes[i].stage == STAGE::INTERC){height = sorted_nodes[i].cluster_num;}
        dendrogram[i] = dendroLine(left, right, height, sorted_nodes[i].size());
    });
    return dendrogram;
}

// updates SP
template<class T, class PROF> 
void DBHTTMFG<T, PROF>::buildHierarchy(){
    auto info = TreeChainInfo(n);
    auto finder = MatrixNNFinder<T, distComplete<T>>(n, &SP);
    auto flags = parlay::sequence<std::size_t>(n);

    std::size_t round = 0;
    auto vtxs = parlay::delayed_seq<std::size_t>(n, [&](std::size_t i){return i;});
    auto sorted_vtx = parlay::sort(vtxs, [&](std::size_t i, std::size_t j){
        pair<std::size_t, std::size_t> a = make_pair(flatClustering[i], bbMember[i]);
        pair<std::size_t, std::size_t> b = make_pair(flatClustering[j], bbMember[j]);
        return a < b;
    });

    //offsets[v] is the starting point of neighbors of v
    auto offsets = parlay::sequence<std::size_t>(n);
    // if i = offsets[v] for some v, sizes[i] stores the size of the neighbor list
    auto sizes = parlay::sequence<std::size_t>(n);

    // intra bubble
    size_t ind = 0;
    offsets[sorted_vtx[0]] = ind;
    for(std::size_t k=1;k<n;++k){
        vtx i = sorted_vtx[k-1];
        vtx j = sorted_vtx[k];
        pair<std::size_t, std::size_t> a = make_pair(flatClustering[i], bbMember[i]);
        pair<std::size_t, std::size_t> b = make_pair(flatClustering[j], bbMember[j]);
        if(a != b){
            ind = k;
            sizes[offsets[i]] = ind - offsets[i];
        }
        offsets[j] = ind;
    }
    sizes[offsets[sorted_vtx[n-1]]] = n- offsets[sorted_vtx[n-1]];
    
    chain_linkage_matrix(&finder, &info, sorted_vtx, offsets, sizes, flags, STAGE::INTRABB, round);

        //intra cluster
    ind = 0;
    offsets[sorted_vtx[0]] = ind;
    for(std::size_t k=1;k<n;++k){
        vtx i = sorted_vtx[k-1];
        vtx j = sorted_vtx[k];
        std::size_t a = flatClustering[i];
        std::size_t b = flatClustering[j];
        if(a != b){
            ind = k;
            sizes[offsets[i]] = ind - offsets[i];
        }
        offsets[j] = ind;
    }
    sizes[offsets[sorted_vtx[n-1]]] = n- offsets[sorted_vtx[n-1]];
    chain_linkage_matrix(&finder, &info, sorted_vtx, offsets, sizes, flags, STAGE::INTERBB, round);

    std::size_t interClusterID = finder.nodeIdx.load();
    //inter cluster
    parlay::parallel_for(0,n,[&](std::size_t i){
        offsets[i]=0;
    });
    sizes[0] = n;
    chain_linkage_matrix(&finder, &info, sorted_vtx, offsets, sizes, flags, STAGE::INTERC, round);

    //output format
    dendro = formatDendrogram(finder.nodes, n, 0, interClusterID, flatClustering.data(), bbMember.data());

}