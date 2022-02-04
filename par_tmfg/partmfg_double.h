#pragma once
#include "partmfg.h"
#include "profiler.h"
#include "dbht.h"

#ifdef PROFILE
typedef ParTMFG<double, Profiler> ParTMFGD;
typedef DBHTTMFG<double, Profiler> ParDBHTTMFGD;
#else
typedef ParTMFG<double, DummyProfiler> ParTMFGD;
typedef DBHTTMFG<double, DummyProfiler> ParDBHTTMFGD;
#endif