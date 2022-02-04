#pragma once


enum STAGE { INTRABB=1, INTERBB=2, INTERC=3, POINT=4 };

struct Node{
    std::size_t cId;
    std::size_t round = 0;
    std::size_t idx; // node idx
    std::size_t n = 1;
    Node *left = nullptr;
    Node *right = nullptr;
    double height = 0;
    std::size_t cluster_num = 1; //number of clusters in this Node, used for inter-cluster
    STAGE stage = POINT;

    Node(std::size_t  t_cid, std::size_t t_round, std::size_t t_idx, Node *t_left, Node *t_right, double _height):
        cId(t_cid),
        round(t_round),
        idx(t_idx),
        left(t_left),
        right(t_right),
        height(_height){
            n = t_left->n + t_right->n;
        }
    
    Node(std::size_t  t_cid):
        cId(t_cid),
        idx(t_cid){
        }
    Node(){}

    inline bool isLeaf() const {return n == 1;}
    inline double getHeight() const {return height;}
    inline std::size_t getRound() const {return round;}
    inline std::size_t getIdx() const {return idx;}
    inline std::size_t size() const {return n;}
};

struct nodeComparator{
  double eps;
  nodeComparator(double _eps):eps(_eps) {
  }

  bool operator() (const Node i, const Node j) const {
      if(abs(i.getHeight() - j.getHeight()) <= eps) return i.getRound() < j.getRound();
     return i.getHeight() < j.getHeight();
  }
  // bool operator() (const Node i, const Node j) const {
  //   if(i.getRound() != j.getRound()) return i.getRound() < j.getRound();
  //   // if(abs(i.getHeight() - j.getHeight()) <= eps) return false;
  //   return i.getHeight() < j.getHeight();
  // }
};


struct EDGE{//nodes unordered
    volatile std::size_t first;
    volatile std::size_t second;
    volatile double w;//weight

    EDGE(std::size_t t_u, std::size_t t_v, double t_w):first(t_u), second(t_v), w(t_w){}
    EDGE():first(NO_NEIGH), second(NO_NEIGH), w(numeric_limits<double>::max()){}

    inline void print(){
        std::cout << "(" << first << ", " <<  second << "):" << w << std::endl;
    }

    inline double getW() const {return w;}
    inline pair<std::size_t, std::size_t> getE() const {return make_pair(first,second);}
    inline void update(std::size_t t_u, std::size_t t_v, double t_w){
      first = t_u;
      second = t_v;
      w = t_w;
    }
};

struct edgeComparator2{
    double eps = 1e-20;
    edgeComparator2(){}
    edgeComparator2(double _eps):eps(_eps){}
    bool operator () (EDGE i, EDGE j) {
      // if(abs(i.getW() - j.getW()) < i.getW()*std::numeric_limits<double>::epsilon()) return i.second < j.second;
      if(abs(i.getW() - j.getW()) <= eps) return i.second < j.second;
      return i.getW() < j.getW();
      }
};


enum Method { WARD, COMP, AVG, AVGSQ };


// same as distComplete1, uses array instead of cache
template<class T>
struct distComplete {
  Method method = COMP;

  inline static void printName(){
    cout << "distComplete" << endl;
  }
  inline double updateDistO(double d1, double d2, double nql, double nqr, double nr, double dij){
    return max(d1,d2);
  }

  inline double updateDistN(double d1, double d2, double d3, double d4, 
                       double nql, double nqr, double nrl, double nrr,
                       double dij, double dklr){
    return max(max(max(d1,d2), d3), d4);
  }
};

  struct dendroLine{
    std::size_t id1;
    std::size_t id2;
    double height;
    std::size_t size;
    dendroLine(std::size_t _id1, std::size_t _id2, double _height, std::size_t _size):id1(_id1), id2(_id2), height(_height), size(_size){}

    void print(ofstream &file_obj){
        file_obj << id1 << " " << id2 << " " << std::setprecision(20) << height << " " << size << endl; 
    }

    void print(){
        cout << id1 << " " << id2 << " " << height << " " << size << endl; 
    }
  };
