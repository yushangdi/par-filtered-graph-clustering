#include "partmfg.cpp"
#include "profiler.h"
#include "dbht_direction.cpp"
#include "dbht_bbtree.cpp"
#include "dbht_bfs.cpp"
#include "dbht_assign.cpp"
#include "dbht_hierarchy.cpp"
#include "dbht_outputs.cpp"


#ifdef PROFILE
template struct ParTMFG<double, Profiler>;
template struct DBHTTMFG<double, Profiler>;
#else
template struct ParTMFG<double, DummyProfiler>;
template struct DBHTTMFG<double, DummyProfiler>;
#endif

