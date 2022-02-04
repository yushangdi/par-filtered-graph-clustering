#pragma once

#include "linkage_types.h"
#include "chain.h"
#include "../common/common.h"

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"

using namespace std;

namespace research_graph {
namespace in_memory {

template<class T, class distF>
class MatrixNNFinder {
  public:
  typedef Node nodeT;

  std::size_t n;// number of points
  std::size_t C;// number of clusters
  constexpr static double eps=0;

  SymMatrix<T> *matrix; // distance to clusters, store two copies of each distance
  distF distComputer;
  edgeComparator2 EC2;

  UnionFind::ParUF<std::size_t> *uf;
  parlay::sequence<std::size_t> activeClusters; // ids of connected components
  parlay::sequence<std::size_t> activeClustersCopy; // ids of connected components
  parlay::sequence<bool> flag;//used in updateActivateClusters, flag[i] is activeClusters[i]
  parlay::sequence<EDGE> edges; //edges[i] stores the min neigh of cluster i
  
  parlay::sequence<std::size_t> rootIdx;
  parlay::sequence<Node> nodes; // preallocated space to store tree nodes
  atomic<std::size_t> nodeIdx; // the index of next node to use for cluster trees

  


  MatrixNNFinder(std::size_t t_n, SymMatrix<T> *t_matrix): 
    n(t_n), C(t_n), matrix(t_matrix){
    distComputer = distF();
    EC2 = edgeComparator2(eps);
    nodeIdx.store(n); // have used n nodes
    uf = new UnionFind::ParUF<std::size_t>(n);

    flag = parlay::sequence<bool>(n); 
    edges = parlay::sequence<EDGE>(n);
    activeClusters = parlay::sequence<std::size_t>(n);
    activeClustersCopy = parlay::sequence<std::size_t>(n);
    rootIdx = parlay::sequence<std::size_t>(n);
    nodes = parlay::sequence<Node>(2*n);
    
    parlay::parallel_for(0,n,[&](std::size_t i){
        nodes[i] = Node(i);
        rootIdx[i] = i;
        edges[i] = EDGE(NO_NEIGH,NO_NEIGH,numeric_limits<double>::max());
        activeClusters[i] = i;
    }) ;
  }

  inline std::size_t cid(nodeT* node){ return node->cId;}
  inline std::size_t idx(nodeT* node){ return node->idx;}
  inline nodeT *getNode(std::size_t cid){return nodes.data()+rootIdx[cid];}
  inline std::size_t idx(std::size_t cid){return idx(getNode(cid));}
  inline std::size_t leftIdx(std::size_t cid){return getNode(cid)->left->idx;}
  inline std::size_t rightIdx(std::size_t cid){return getNode(cid)->right->idx;}
  inline std::size_t cid(std::size_t idx){return cid(nodes[idx]);}
  inline bool isActive(std::size_t cid){return uf->find(cid)==cid;}

  // i, j are cluster ids
  // find distance in matrix
  double getDist(std::size_t i,  std::size_t j){
    return matrix->get(i, j);
  }

  double getDist(nodeT *inode,  nodeT *jnode){
    return getDist(cid(inode), cid(jnode));
  }

  //newc is a newly merged cluster
  // newc is new, rid is old
  inline double getNewDistO(std::size_t newc, std::size_t rid){
    nodeT* ql = getNode(newc)->left;
    std::size_t nql = ql->n;
    nodeT* qr = getNode(newc)->right;
    std::size_t nqr = qr->n;
    double dij = getNode(newc)->getHeight();
    nodeT* rroot = getNode(rid);
    std::size_t nr = rroot->n;

    double d1,d2;
    d1 = getDist(ql,rroot);
    d2 = getDist(qr,rroot);

    return distComputer.updateDistO(d1, d2, nql, nqr, nr, dij);
  }


  //newc is a newly merged cluster
  // newc is new, rid is merged
  inline double getNewDistN(std::size_t newc, std::size_t rid){
    nodeT* ql = getNode(newc)->left;
    std::size_t nql = ql->n;
    nodeT* qr = getNode(newc)->right;
    std::size_t nqr = qr->n;
    double dij = getNode(newc)->getHeight();

    nodeT* rl = getNode(rid)->left;
    std::size_t nrl = rl->n;
    nodeT* rr = getNode(rid)->right;
    std::size_t nrr = rr->n;
    double dklr = getNode(rid)->getHeight();

    double d1,d2, d3, d4;
    d1 = getDist(ql,rl);
    d2 = getDist(ql,rr);
    d3 = getDist(qr,rl);
    d4 = getDist(qr,rr);

    return distComputer.updateDistN(d1, d2, d3, d4, nql, nqr, nrl, nrr, dij, dklr);
  }

  struct neighborComparator{
    std::size_t cid;
    SymMatrix<T> *matrix;
    neighborComparator(std::size_t _cid, SymMatrix<T> *_matrix):cid(_cid), matrix(_matrix){}
    pair<double, std::size_t> get(std::size_t i) {
        if(i==cid) return make_pair(numeric_limits<double>::max(), i); //used for max index, so reverse comparason
        return make_pair(matrix->get(i, cid), i);
    }

    bool operator() (std::size_t a, std::size_t b) {
        pair<double, std::size_t> i  = get(a);
        pair<double, std::size_t> j  = get(b);
        double dist1 = i.first;
        double dist2 = j.first;
        if(abs(dist1 - dist2) <= eps) return i.second < j.second;
        return dist1 < dist2;
    }

    };

  // store the closest nn among neighbors in edges
  // neighbors is a list of active cluster ids
  inline void getNN(std::size_t cid, parlay::slice<std::size_t*, std::size_t*> &neighbors){
    std::size_t* nn = parlay::min_element(neighbors, neighborComparator(cid, matrix));
    double minD = matrix->get(cid, *nn);
    edges[cid].update(cid, *nn, minD);
  }

    // merge two tree nodes with cID u and v into a new tree node with cId newc
    // u, v are cluster ids
  inline void merge(std::size_t u, std::size_t v, std::size_t newc, std::size_t round, double height){
    std::size_t rootNodeIdx = nodeIdx.fetch_add(1);
    nodes[rootNodeIdx] = nodeT(newc, round, rootNodeIdx, getNode(u), getNode(v), height);
    rootIdx[newc] = rootNodeIdx;
  }

  inline void updateDist(std::size_t newc, std::size_t round){
    // std::size_t idx1 = leftIdx(newc);
    // std::size_t idx2 = rightIdx(newc);
    std::size_t cid1 = getNode(newc)->left->cId;
    std::size_t cid2 = getNode(newc)->right->cId;

    // std::size_t oldc = newc == cid1? cid2:cid1;
    parlay::parallel_for(0, C, [&](std::size_t i){
        std::size_t updc = activeClusters[i];
        if(updc != cid1 && updc != cid2 && uf->find(updc) == updc){
          if(getNode(updc)->round == round){ /* if updc is also a merged cluster */
            if(updc > newc){ /* the smaller id cluster will update */
                double dist = getNewDistN(newc, updc); 
                matrix->update(updc, newc, dist);
            }
          }else{
              double dist = getNewDistO(newc, updc);
              matrix->update(updc, newc, dist);
          }

        }//end of (updc != cid1 && updc != cid2)
    });//end parallel for
  }

  // find new activeClusters array based on uf, update C
  inline void updateActiveClusters(std::size_t round){
      parlay::parallel_for(0,C,[&](std::size_t i){
        std::size_t cid  = activeClusters[i];
        flag[i] = isActive(cid);
        activeClustersCopy[i] = activeClusters[i];
      });
    C = parlay::pack_into(make_slice(activeClustersCopy).cut(0,C), flag, activeClusters);
  }
}; // finder end


vector<dendroLine> formatDendrogram(parlay::sequence<Node> &nodes, std::size_t n, double eps){
    auto sorted_nodes = parlay::sort(parlay::make_slice(nodes).cut(n,2*n-1), nodeComparator(eps));

    auto map = parlay::sequence<std::size_t>(n);
    parlay::parallel_for(0, n-1, [&](std::size_t i){
        map[sorted_nodes[i].getIdx() - n] = i+n;
    });
    vector<dendroLine> dendrogram = vector<dendroLine>(n-1);//(dendroLine *)malloc(sizeof(dendroLine)*(n-1));
    parlay::parallel_for(0, n-1, [&](std::size_t i){
        std::size_t left = sorted_nodes[i].left->getIdx();
        std::size_t right = sorted_nodes[i].right->getIdx();

        left = left < n ? left : map[left-n];
        right = right < n ? right : map[right-n];

        if(left > right) swap(left, right);


        dendrogram[i] = dendroLine(left, right, sorted_nodes[i].getHeight(),sorted_nodes[i].size());
    });
    return dendrogram;
}


template<class TF>
inline void chain_find_nn(std::size_t chainNum, TF *finder, TreeChainInfo *info){
  auto c_neigh = make_slice(finder->activeClusters).cut(0, finder->C);
  // update edges array
  parlay::parallel_for(0, chainNum, [&](std::size_t i){
    std::size_t cid = info->terminal_nodes[i];
    // finder->edges[cid] = EDGE(cid,NO_NEIGH,numeric_limits<double>::max());
    finder->getNN(cid, c_neigh);
  });
  // update chain
  parlay::parallel_for(0, chainNum, [&](std::size_t i){
      std::size_t cid = info->terminal_nodes[i];
      info->updateChain(cid, finder->edges[cid].second, finder->edges[cid].getW());
  });
}


template<class TF>
inline void link_terminal_nodes(TF *finder, TreeChainInfo *info, std::size_t round, parlay::sequence<std::size_t> &flags){
  std::size_t chainNum = info->chainNum;
  auto uf = finder->uf;
  EDGE *edges = finder->edges.data();

  parlay::parallel_for(0, chainNum, [&](std::size_t i){
    std::size_t cid = info->terminal_nodes[i];
    if(edges[cid].second != NO_NEIGH){ // do not link if there no neighbor in this round
        std::size_t nn = edges[cid].second; 
        std::size_t nn_of_nn =info->getNN(nn);

        // avoid merging twice
        // update chain info and merge trees
        if ((!(info->isTerminal(nn) && cid <= nn)) && (nn_of_nn == cid)){// if they are RNN 
        std::size_t newc = uf->link(cid, nn, edges[cid].getW());
        info->invalidate(nn, NO_NEIGH);
        info->invalidate(cid, NO_NEIGH);
        finder->merge(cid, nn, newc, round, edges[cid].getW());
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


template<class T, class distT>
vector<dendroLine> chain_linkage_matrix(SymMatrix<T>* M){
  std::size_t n = M->n;
  auto info = TreeChainInfo(n);
  auto finder = MatrixNNFinder<T, distT>(n, M);
  std::size_t chainNum = info.chainNum;
  int round = 0;
  auto flags = parlay::sequence<std::size_t>(n);
  while(finder.C > 1){
    round ++;
#ifdef DEBUG
    if(round >= n){
        cout << "too many rounds" << endl;
        exit(1);
    }
#endif
    chain_find_nn(chainNum, &finder, &info);
    link_terminal_nodes<MatrixNNFinder<T, distT>>(&finder, &info, round, flags);
    // get ready for next round
    finder.updateActiveClusters(round);
    info.next(&finder);
    chainNum = info.chainNum;
  }
  return formatDendrogram(finder.nodes, n, 0);
}


}  // namespace in_memory
}  // namespace research_graph
