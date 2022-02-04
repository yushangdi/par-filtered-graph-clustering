#pragma once

#include "unionfind.h"
#include "linkage_types.h"



template<class T, class distF>
class MatrixNNFinder {
  public:
  typedef Node nodeT;

  std::size_t n;// number of points
  std::size_t C;// number of clusters
  constexpr static double eps=0;

  SymM<T> *matrix; // distance to clusters, store two copies of each distance
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

  


  MatrixNNFinder(std::size_t t_n, SymM<T> *t_matrix): 
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
    SymM<T> *matrix;
    neighborComparator(std::size_t _cid, SymM<T> *_matrix):cid(_cid), matrix(_matrix){}
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
  inline void getNN(std::size_t cid, parlay::slice<std::size_t*, std::size_t*> &neighbors){//parlay::slice<std::size_t>
    std::size_t* nn = parlay::min_element(neighbors, neighborComparator(cid, matrix));
    double minD = matrix->get(cid, *nn);
    edges[cid].update(cid, *nn, minD);
  }

    // merge two tree nodes with cID u and v into a new tree node with cId newc
    // u, v are cluster ids
  inline std::size_t merge(std::size_t u, std::size_t v, std::size_t newc, std::size_t round, double height){
    std::size_t rootNodeIdx = nodeIdx.fetch_add(1);
    nodes[rootNodeIdx] = nodeT(newc, round, rootNodeIdx, getNode(u), getNode(v), height);
    rootIdx[newc] = rootNodeIdx;
    return rootNodeIdx;
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
