#pragma once

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "../common/atomics.h"
#include "linkage_types.h"

using namespace std;

namespace research_graph {
namespace in_memory {

// *************************************************************
//   NN Chain
// *************************************************************
struct TreeChainInfo{
  parlay::sequence<std::size_t> terminal_nodes;// must be cluster id
  parlay::sequence<std::size_t> chain;//chain[i] is cluster i's nn, NO_NEIGH for unknown or invalid// -1 for unknown, -2 for invalid
  parlay::sequence<bool> is_terminal;
  parlay::sequence<bool> flag;
  std::size_t chainNum;

  TreeChainInfo(std::size_t n){
    terminal_nodes = parlay::sequence<std::size_t>(n);
    parlay::parallel_for(0,n,[&](std::size_t i){terminal_nodes[i] = i;});
    chain = parlay::sequence<std::size_t>(n, NO_NEIGH);
    is_terminal = parlay::sequence<bool>(n, true);
    flag = parlay::sequence<bool>(n);
    chainNum = n;
  }

  inline void updateChain(std::size_t cid, std::size_t nn, double w){
    chain[cid] = nn;
  }
  // then in checking, only check for -1
  inline void invalidate(std::size_t cid, std::size_t code){
    chain[cid] = code;
  }
  inline std::size_t getNN(std::size_t cid){return chain[cid];}
  inline bool isTerminal(std::size_t cid){return is_terminal[cid];}

  // update findNN, terminal_nodes and chainNum
  template<class F>
  inline void next(F *finder){
    parlay::parallel_for(0, finder->n, [&](std::size_t i){
      is_terminal[i] = false;
    });
    std::size_t C = finder->C;
    parlay::parallel_for(0, C, [&](std::size_t i){
      std::size_t cid = finder->activeClusters[i];
      flag[i] = getNN(cid) == NO_NEIGH || getNN(getNN(cid)) == NO_NEIGH ;// only merged clusters have negative neighbor in chain ok because -2 won't be in active clusters
      is_terminal[cid] = flag[i];
    });
    chainNum = parlay::pack_into(make_slice(finder->activeClusters).cut(0,C), flag, terminal_nodes);
  }

};

}  // namespace in_memory
}  // namespace research_graph