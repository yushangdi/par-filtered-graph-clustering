#include <string>
#include <vector>

#include "linkage.h"
#include "../common/common.h"

using namespace std;

int main(int argc, char *argv[]) {
    char* filename = argv[1];
    std::size_t n = atoi(argv[2]);
    std::string output = argv[3];
    std::string method = argv[4];
    int round = atoi(argv[5]);

    cout << "num workers: " << parlay::num_workers() << endl;

    research_graph::in_memory::SymMatrix<double> W = research_graph::in_memory::IO::readSymMatrixFromFile<double>(filename, n); 

    parlay::parallel_for(0, W.distmat_size, [&](size_t i){
        W.distmat[i] = sqrt(2*(1-W.distmat[i]));
    });


    timer t1;t1.start();
    vector<research_graph::in_memory::dendroLine> dendro;

    for(int r=0;r<round;++r){
    cout << "=========" << endl;
    if(method == "comp"){
        using distT = research_graph::in_memory::distComplete<double>;
        dendro = research_graph::in_memory::chain_linkage_matrix<double, distT>(&W);
    }else if(method == "avg"){
        using distT = research_graph::in_memory::distAverage<double>;
        dendro = research_graph::in_memory::chain_linkage_matrix<double, distT>(&W);
    }else{
        cout << "invalid method " << method << endl;
        exit(1);
    }
    cout << "hac_time: " << t1.next() << endl;
    }
    

    ofstream file_obj;
    file_obj.open(output); 
    for(size_t i=0;i<n-1;i++){
        dendro[i].print(file_obj);
    }
    file_obj.close();
    cout << "done" << endl;

}
