#ifndef EFANNA2E_INDEX_NSG_H
#define EFANNA2E_INDEX_NSG_H

#include "util.h"
#include "parameters.h"
#include "neighbor.h"
#include "index.h"
#include <cassert>
#include <unordered_map>
#include <string>
#include <sstream>
#include <boost/dynamic_bitset.hpp>
#include <stack>

    typedef unsigned short int vl_type;

    class VisitedList {
    public:
        vl_type curV;
        vl_type *mass;
        unsigned int numelements;

        VisitedList(int numelements1) {
            curV = -1;
            numelements = numelements1;
            mass = new vl_type[numelements];
        }

        void reset() {
            curV++;
            if (curV == 0) {
                memset(mass, 0, sizeof(vl_type) * numelements);
                curV++;
            }
        };

        ~VisitedList() { delete[] mass; }
    };
///////////////////////////////////////////////////////////
//
// Class for multi-threaded pool-management of VisitedLists
//
/////////////////////////////////////////////////////////

    class VisitedListPool {
        std::deque<VisitedList *> pool;
        std::mutex poolguard;
        int numelements;

    public:
        VisitedListPool(int initmaxpools, int numelements1) {
            numelements = numelements1;
            for (int i = 0; i < initmaxpools; i++)
                pool.push_front(new VisitedList(numelements));
        }

        VisitedList *getFreeVisitedList() {
            VisitedList *rez;
            {
                std::unique_lock <std::mutex> lock(poolguard);
                if (pool.size() > 0) {
                    rez = pool.front();
                    pool.pop_front();
                } else {
                    rez = new VisitedList(numelements);
                }
            }
            rez->reset();
            return rez;
        };

        void releaseVisitedList(VisitedList *vl) {
            std::unique_lock <std::mutex> lock(poolguard);
            pool.push_front(vl);
        };

        ~VisitedListPool() {
            while (pool.size()) {
                VisitedList *rez = pool.front();
                pool.pop_front();
                delete rez;
            }
        };
    };


struct elem{
	unsigned int id;
	float dist;
	bool flag;
};

namespace efanna2e {

class IndexNSG : public Index {
 public:
  explicit IndexNSG();
  explicit IndexNSG(const size_t dimension, const size_t n, Metric m, Index *initializer, const bool flag, float*** book, const int l, float*** quanti, int level);


  virtual ~IndexNSG();

  virtual void Save(const char *filename)override;
  virtual void Load(const char *filename)override;
  virtual void Init(size_t dimension, size_t n, Metric m, Index *initializer, 
        const bool flag, float*** book, int l, float*** quanti);

  virtual void Build(size_t n, float *data, const Parameters &parameters);
  virtual void Build2(size_t n, unsigned char *data, const Parameters &parameters);

  virtual void Search(
      const float *query,
      float *x,
      size_t k,
      const Parameters &parameters,
      unsigned *indices) ;
  void SearchWithOptGraph(
      const float *query,
      size_t K,
      const Parameters &parameters,
      unsigned *indices,
      unsigned *init_point,
	  VisitedListPool *visited_list_pool_);
  void SearchWithquanGraph(float** book_, int X, unsigned* points, elem* init_point, 
        unsigned int* trans, const Parameters &parameters, char* fflag, int level, VisitedListPool *visited_list_pool_); 
  void SearchWithquanGraph2(float** book_, int X, elem* points, unsigned* init_point, 
        unsigned int* trans, const Parameters &parameters, VisitedListPool *visited_list_pool_); 			   
		
  void SearchWithquanGraph3(float** book_, int X, elem* points, elem* init_point, 
        unsigned int* trans, const Parameters &parameters, char* fflag, int level, VisitedListPool *visited_list_pool_);
	
  void SearchWithsingleGraph(float** book_, int X, unsigned* points, unsigned* enter_obj, 
        unsigned int* trans, const Parameters &parameters, VisitedListPool *visited_list_pool_);

	
  void OptimizeGraph(float* data);
  void OptimizeQuan(unsigned char *quan_, int level, int num, int max_level);
  char** quan_graph;
  
  protected:
    typedef std::vector<std::vector<unsigned > > CompactGraph;
    typedef std::vector<SimpleNeighbors > LockGraph;
    typedef std::vector<nhood> KNNGraph;

    CompactGraph final_graph_;

    Index *initializer_;
    void init_graph(const Parameters &parameters);
	void init_graph2(const Parameters &parameters);
    void get_neighbors(
        const float *query,
        const Parameters &parameter,
        std::vector<Neighbor> &retset,
        std::vector<Neighbor> &fullset);
    void get_neighbors(
        const float *query,
        const Parameters &parameter,
        boost::dynamic_bitset<>& flags,
        std::vector<Neighbor> &retset,
        std::vector<Neighbor> &fullset);
    //void add_cnn(unsigned des, Neighbor p, unsigned range, LockGraph& cut_graph_);
    void InterInsert(unsigned n, unsigned range, std::vector<std::mutex>& locks, SimpleNeighbor* cut_graph_);
    void sync_prune(unsigned q, std::vector<Neighbor>& pool, const Parameters &parameter, boost::dynamic_bitset<>& flags, SimpleNeighbor* cut_graph_);
    void Link(const Parameters &parameters, SimpleNeighbor* cut_graph_);
    void Load_nn_graph(const char *filename);
    void tree_grow(const Parameters &parameter);
    void DFS(boost::dynamic_bitset<> &flag, unsigned root, unsigned &cnt);
    void findroot(boost::dynamic_bitset<> &flag, unsigned &root, const Parameters &parameter);


  private:
    unsigned width;
    unsigned ep_;
    std::vector<std::mutex> locks;
    char* opt_graph_;
    size_t node_size;
    size_t data_len;
    size_t neighbor_len;
    KNNGraph nnd_graph;
	//-----------added-------------
	//size_t quan_len[max_level];
	//size_t quan_size[max_level];
		
	//-----------------------------
};
}

#endif //EFANNA2E_INDEX_NSG_H
