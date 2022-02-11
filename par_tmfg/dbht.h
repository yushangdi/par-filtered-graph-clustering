#pragma once

#include "macro.h"
#include <tuple>
#include <limits>
#include <atomic>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/delayed_sequence.h"

#include "gettime.h"
#include "atomics.h"
#include "utilities.h"
#include "symmatrix.h"
#include "linkage_types.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace std;


struct tbbTree{ //TMFG bubble Tree
    std::size_t parent = DBHT_NONE;
    std::size_t children[3] = {DBHT_NONE, DBHT_NONE, DBHT_NONE};
    atomic<int> k=0;

    void init(){
        parent = DBHT_NONE;
        for(int i=0;i<3;++i){ children[i]=DBHT_NONE;}
        k.store(0);
    }

    void addChild(std::size_t i){
#ifdef DEBUG
        if(k>=3){cout << "adding a fourth child" << endl;}
#endif
        children[k] = i; k++;
    }

    void addChildPar(std::size_t i){
        int kk = k.fetch_add(1);
#ifdef DEBUG
        if(kk>=3){cout << "adding a fourth child" << endl;}
#endif
        children[kk] = i;
    }

    // inline bool isEmpty(int kk){return children[kk]==DBHT_NONE;}
    inline bool isEmpty(int kk){return kk>=k;}
};


template<class T, class PROF> 
struct DBHTTMFG{
    cliqueT* cliques;
    triT* triangles; 
    std::size_t n;
    SymM<T>* W;
    SymM<T>* D;
    tuple<vtx,vtx,T> *P;
    PROF *pf;   
    std::size_t *t2c; // triangle to cliqiue
    tbbTree *tree;
    T *degrees;
    bool *directions; // True is pointing to parent, False is pointing from parent
    parlay::sequence<bool> converging;
    std::size_t root; // index of root clique in cliques

    SymM<T> SP;// table of shortest path
    parlay::sequence<bool> v2cl;
    // from vertex to Converging bubble. 
    // false if the vertex cannot reach the conveging bubble, true if can reach
    parlay::sequence<bool> b2cl; //from bubble to converging bubble

    std::size_t nb; // number of bubbles
    std::size_t nc; // numbner of converging bubbles
    parlay::sequence<std::size_t> cvg_bubbles;
    parlay::sequence<std::size_t> map;//map from bubble/clique id to location in cvg_bubbles

    parlay::sequence<std::size_t> flatClustering; //cluster membership of vertices
    parlay::sequence<std::size_t> bbMember; //bubble membership of vertices
    parlay::sequence<bool> assigned;
    dendroLine* dendro;

    // cliques should have been inited
    DBHTTMFG(cliqueT* _cliques, triT* _triangles, std::size_t _n, SymM<T> *_W, tuple<vtx,vtx,T> *_P, SymM<T> *_D, PROF *_profiler):
    cliques(_cliques), triangles(_triangles), n(_n), W(_W), D(_D), P(_P), pf(_profiler){
        nb = n-3;
        root = 0;
        t2c = (std::size_t *)malloc(3*_n * sizeof(std::size_t)); //init others?
        t2c[0] = 0;t2c[1] = 0;t2c[2] = 0;t2c[3] = 0;

        tree = (tbbTree *)malloc(_n * sizeof(tbbTree));
        degrees = (double *)malloc(n * sizeof(double));
        directions = (bool *)malloc((n-3) * sizeof(bool)); //init?
        converging = parlay::sequence<bool>(n, true);

        parlay::parallel_for(0, _n, [&](std::size_t i) {
            tree[i].init();
            degrees[i] = 0;
        });

        SP.init(n);
        map = parlay::sequence<std::size_t>(n);
        flatClustering = parlay::sequence<std::size_t>(n, 0);
        bbMember = parlay::sequence<std::size_t>(n);

        initDegrees();
    }

    ~DBHTTMFG(){
        free(tree);
        free(t2c);
        free(degrees);
        free(directions);
        free(dendro);
        // free(v2c);
        // free(converging);
    }

    inline T getW(vtx i, vtx j){
        return W->get(i,j);
    }

    inline T getD(vtx i, vtx j){
        if(D==nullptr){
            if(i==j) return 0;
            T w = getW(i,j);
            return sqrt(2*(1-w));
        }
        return D->get(i,j);
    }
    

    void initDegrees(){
        vtx t1, t2, t3,t4;
        tie(t1,t2,t3,t4) = cliques[0];
        degrees[t1] += getW(t1, t2) + getW(t1, t3) + getW(t1, t4);
        degrees[t2] += getW(t2, t1) + getW(t2, t3) + getW(t2, t4);
        degrees[t3] += getW(t3, t1) + getW(t3, t2) + getW(t3, t4);
        degrees[t4] += getW(t4, t1) + getW(t4, t2) + getW(t4, t3);
    }
    
    void APSP(){
        // timer t1;t1.start();
        typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
            boost::no_property, boost::property < boost::edge_weight_t, T > > graph_t;
        typedef typename boost::graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
        typedef std::pair<int, int> Edge;

        // //prepare adj matrix, inlcude zero weight edges
        // const std::size_t num_nodes = n;
        // std::size_t num_arcs = 3*n-6;
        // vector<Edge> edge_array = vector<Edge>(num_arcs);
        // vector<T> weights = vector<T>(num_arcs);
        // for(size_t i=0;i<num_arcs;i++){ //TODO: parallel_for
        //     vtx a = get<0>(P[i]); vtx b = get<1>(P[i]);
        //     edge_array[i] = Edge(a, b);
        //     weights[i] = getD(a, b);
        // }
        // graph_t g(edge_array.begin(), edge_array.end(), weights.begin(), num_nodes);

        //prepare adj matrix, not inlcude zero weight edges
        const std::size_t num_nodes = n;
        std::size_t num_arcs = 3*n-6;
        auto inds = parlay::delayed_seq<bool>(num_arcs, [&](std::size_t i){
            return getD(get<0>(P[i]), get<1>(P[i])) > 0;
        });
        auto pos_inds = parlay::pack_index(inds);
        num_arcs = pos_inds.size();

        vector<Edge> edge_array = vector<Edge>(num_arcs);
        vector<T> weights = vector<T>(num_arcs);
        parlay::parallel_for(0, num_arcs, [&](size_t i){ //for(size_t i=0;i<num_arcs;i++){
            vtx a = get<0>(P[pos_inds[i]]); vtx b = get<1>(P[pos_inds[i]]);
            edge_array[i] = Edge(a, b);
            weights[i] = getD(a, b);
        });
        graph_t g(edge_array.begin(), edge_array.end(), weights.begin(), num_nodes);

        // cout << "build: " << t1.next() << endl;

        // compute to vertex v
        parlay::parallel_for(0, n, [&](vtx v){
            std::vector<vertex_descriptor> p(num_vertices(g));
            std::vector<T> d(num_vertices(g));
            vertex_descriptor s = vertex(v, g);
            boost::dijkstra_shortest_paths(g, s, boost::predecessor_map(&p[0]).distance_map(&d[0]));

            parlay::parallel_for(v, n, [&](vtx u){
                SP.update(v, u, d[u]);
            });
        });
    }

    // inserting v to triangle [v1, v2, v3]
    void updateDegrees(vtx t1, vtx t2, vtx t3, vtx v);

    // t1 is the triangle replaced, t2, t3 are the two other new triangles' indices in triangles
    // newc is the index of the new clique in cliques
    // v is the new vertex inserted
    void insertOne(std::size_t t1, std::size_t t2, std::size_t t3, std::size_t newc, vtx v);


    //////// direction //////////////////
    struct partialDeg{
        T vals[3]={0};  
        vtx inds[3];

        void inc(vtx i, T val){
            for(int k=0;k<3;++k){
                if(inds[k] == i){
                    pbbs::write_add(&vals[k], val);
                    break;
                }
            }
        }

        void inc(partialDeg &cresult){
            for(int k=0;k<3;++k){
                inc(cresult.inds[k], cresult.vals[k]);
            }
        }
    };

    void computeDirection();

    // return the degree sum inside the triangle (not inclusive) of the three vertex
    // id: the clique id of the current clique, cannot be the root
    // pid: clique id  of the parent
    partialDeg computeDirectionHelper(std::size_t id, std::size_t pid);

    template<int k>
    bool cChelper(std::size_t id1, std::size_t id2){
        return (get<k>(cliques[id1]) == get<0>(cliques[id2]) || 
                get<k>(cliques[id1]) == get<1>(cliques[id2]) || 
                get<k>(cliques[id1]) == get<2>(cliques[id2]) || 
                get<k>(cliques[id1]) == get<3>(cliques[id2]) );
    }
    // compare two cliques with id1 and id2 in cliques
    // prereq: id1 and id2 are parent-child relationship
    //          so three of the elements are the same, and one is different
    // the same ones are at the beginning and the different one (from id1) is in the end
    cliqueT compareCliques(std::size_t id1, std::size_t id2);



    //////// bubble clustering //////////////////

    //bfs on the bubble tree, starting from the bb bubble
    // need to be called after direction is computed
    void bfs(std::size_t bb, bool *visited, std::size_t nc, std::size_t* map);

    // need to be called after direction is computed
    // fill the map and c2cl arrays
    void nonDiscreteClustering();


    struct AT{ 
        T first = NO_ATTACH_MIN; //attachment to the bubble/clique
        std::size_t second=SIZE_T_MAX; //id of the bubble/clique
        AT(){}
        AT(T i, std::size_t j):first(i),second(j){}
    };
    struct cAT{
        bool operator () (const AT i, const AT j){
            if(i.first == j.first) return i.second < j.second;
            return i.first < j.first;
        }
    };

    // assign the vertices to converging bubble

    // compute the attachement from u to bubble i
    // do not consider bubble size, because all bubbles have size 4
    //TODO: need the size in second assignment phase?
    // pereq: u must be a vertex in bubble i
    inline T computeChi(vtx u, std::size_t i){
        T chi = 0;

        vtx v = get<0>(cliques[i]);
        if(v!= u) chi+=getW(u, v);
        v = get<1>(cliques[i]);
        if(v!= u) chi+=getW(u, v);
        v = get<2>(cliques[i]);
        if(v!= u) chi+=getW(u, v);
        v = get<3>(cliques[i]);
        if(v!= u) chi+=getW(u, v);
        // do not divide by the bubble size, because all bubbles have size 4
        return chi;
    }
    void assignToConvergingBubble(); // need to test
    void assignToBubble(); // need to test



    //////////////// build hierarchy ///////////////////////////

    void buildHierarchy();


    //////////////// output ///////////////////////////
    void outputDirections();
    void outputDefaultClustering(string filename, size_t _n=0);
    void outputBubbleAssign(string filename, size_t _n=0);
    void outputDendro(string filename, size_t _n=0);

#ifdef DEBUG
    void checkDegree(){
        for(int i=0; i<n;++i){
            if(degrees[i]<0 || degrees[i] > 8){
                cout << "check degree " << i << endl;
            }
        }
    }
#endif
};
