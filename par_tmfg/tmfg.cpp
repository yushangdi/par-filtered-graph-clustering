// #include "dbht.h"
#include "partmfg_double.h"
#include "profiler.h"
#include "IO.h"

void runDBHT(SymM<double> *W, SymM<double> *D, size_t n, size_t THRESHOLD, string method, bool use_heap, string dsname = ""){
cout << "====" << endl;
cout << "threshold: " << THRESHOLD << endl;
cout << "method: " << method << endl;
cout << "use_heap: " << use_heap << endl;
#ifdef PROFILE
    auto pf = Profiler();
#else
    auto pf = DummyProfiler();
#endif
timer t2;t2.start();
ParTMFGD computer = ParTMFGD(W, n, &pf, use_heap); 
timer t;t.start();
computer.init();           
computer.initGainArray();                 pf.setInitTime(t2.next());
auto clusterer = new ParDBHTTMFGD(computer.cliques.data(), computer.triangles.data(), n, computer.W, computer.P.data(), D, &pf);
int round=0;

if(method == "prefix"){ //get best gain by scanning vertex list
    while(computer.hasUninsertedV()){
            size_t round_THRESHOLD = min(THRESHOLD, computer.getTrianglesNum());
            auto insert_list = computer.getBestVertices(round_THRESHOLD);   pf.incVTime(t2.next());
            computer.inertMultiple(insert_list, clusterer);                 pf.incInsertTime(t2.next());
            computer.updateGainArray(insert_list);                          pf.incUpdTime(t2.next());
#ifdef DEBUG
        computer.checkTriangles();
#endif
            round++;
    } //while end
}else if(method == "exact"){ //use exact
    while(computer.hasUninsertedV()){
        computer.insertOne(clusterer);
        round++;
    } //while end
}else if(method == "naive"){ //naive method
    while(computer.hasUninsertedV()){
        auto insert_list = computer.getAllBestVertices(computer.getTrianglesNum()); pf.incVTime(t2.next());
        computer.inertMultiple(insert_list, clusterer);                             pf.incInsertTime(t2.next());
        computer.initGainArray();                                                   pf.incUpdTime(t2.next());
        round++;
    } //while end
}
// compute the total gain here
cout << "tmfg total: "<< t.next() << endl;
cout << "round: " << round << endl;
computer.computeCost();
pf.report();

t.next();
clusterer->APSP();
cout << "APSP total: "<< t.next() << endl;
clusterer->computeDirection();
cout << "direction total: "<< t.next() << endl;
clusterer->nonDiscreteClustering();
cout << "non-discrete total: "<< t.next() << endl;
clusterer->assignToConvergingBubble(); // need to test
cout << "discrete total: "<< t.next() << endl;
cout << "num cluster: "<< clusterer->nc << endl;
clusterer->assignToBubble(); // need to test
cout << "bubble total: "<< t.next() << endl;
clusterer->buildHierarchy();
cout << "hierarchy total: "<< t.next() << endl;


if(method == "exact" || method == "naive"){
    computer.outputP(dsname + "-" + method + "-P-1");
    clusterer->outputDendro(dsname + "-" + method + "-Z-1" );
}else{
    computer.outputP(dsname + "-" + method + "-P-" + to_string(THRESHOLD) );
    clusterer->outputDendro(dsname + "-" + method + "-Z-" + to_string(THRESHOLD));
}

}

int main(int argc, char *argv[]) {
char* filename = argv[1];
bool use_heap = 0;
string dsname = argv[2];
size_t n = atoi(argv[3]);
char* distance_filename = argv[4]; //if 0, use sqrt(2*(1-w)) by default
string method = argv[5];
size_t THRESHOLD = atoi(argv[6]); // the prefix to insert, only used when method='prefix
int round = atoi(argv[7]);

cout << "workers: " << parlay::num_workers() << endl;
timer t2;t2.start();
SymM<double> W = IO::readSymMatrixFromFile<double>(filename, n); 
W.setDiag(1);
for(int r=0;r<round;++r){
    if(distance_filename[0]!='0'){
        SymM<double> D = IO::readSymMatrixFromFile<double>(distance_filename, n);
        cout << "read: " << t2.next() << endl;
        runDBHT(&W, &D, n, THRESHOLD, method, use_heap, dsname);
    }else{
        cout << "read: " << t2.next() << endl;
        runDBHT(&W, nullptr, n, THRESHOLD, method, use_heap, dsname);
    }
}
}
