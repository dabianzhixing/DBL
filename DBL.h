#include <vector>
#include <tuple>
#include <fstream>
#include <bitset>
#include <set>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/bitset.hpp>
#include <boost/serialization/set.hpp>
#define dl 64 
#define bl 64 
#define query_threadNum 1 
#define construct_threadNum 1

#define _CHECK(expr)                                                \
  if (expr) {                                                           \
  } else {                                                              \
    fprintf(stderr, "CHECK Failed (%s:%d): %s\n",                       \
            __FILE__, __LINE__, #expr);                                 \
    exit(EXIT_FAILURE);                                                 \
  }




namespace boost {
namespace serialization {
class access;
}
}

class DBL {
 public:

  DBL() : V(0) {}


  void construct_index(char *filename1);
  void construct_index(std::istream &ifs1);
  void construct_index_dl(
			std::vector<int> &csr,
			std::vector<int> &prefixsum,
			std::vector<int> &rev_csr,
			std::vector<int> &rev_prefixsum
			);
  void construct_index_bl(
		  	std::vector<int> &csr,
			std::vector<int> &prefixsum,
			std::vector<int> &rev_csr,
			std::vector<int> &rev_prefixsum
			);

  void get_label();
  void insert_edges(std::vector<std::tuple<int,int>> &es);
  void insert_edges(const char *filename, int num);
  void insert_edges(std::istream &ifs, int num);


  void query(
  		std::vector<int> &u,
  		std::vector<int> &v,
  		int sz,
  		bool *result
  		);

  void query_label(
  		std::vector<int> &u,
  		std::vector<int> &v,
  		int sz,
  		bool *sure,
  		bool *result
  		);
  void test_query(
 		std::vector<int> &u,
 		std::vector<int> &v,
 		int sz,
 		bool *test_result,
 		bool *result
 		);

  int V;
  int E;
  int csr_size;
  std::vector<std::set<int> > adj;
  std::vector<std::set<int> > rev_adj;
  std::vector<int> ord;
  std::vector<int> rev_ord;
  std::vector<std::bitset<dl> > dl_in;
  std::vector<std::bitset<dl> > dl_out;
  std::vector<std::bitset<bl> > bl_in;
  std::vector<std::bitset<bl> > bl_out;
  std::vector<int> scc;
  std::vector<int> csr;
  std::vector<int> prefixsum;
  std::vector<int> rev_csr;
  std::vector<int> rev_prefixsum;
  std::vector<int> dl_hash; 
  void get_root_order(std::vector<int> &ord,std::vector<int> &rev_ord);
  void get_root_order_max(std::vector<int> &ord,std::vector<int> &rev_ord);
  void get_root_order_min(std::vector<int> &ord,std::vector<int> &rev_ord);
  void get_root_order_plus(std::vector<int> &ord,std::vector<int> &rev_ord);


  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, unsigned int ver __attribute__((unused))) 
  {
    ar & V;
    ar & E;
    ar & csr_size;
    ar & adj;
    ar & rev_adj;
    ar & ord;
    ar & rev_ord;
    ar & dl_in;
    ar & dl_out;
    ar & bl_in;
    ar & bl_out;
    ar & scc;
    ar & csr;
    ar & prefixsum;
    ar & rev_csr;
    ar & rev_prefixsum;
    ar & dl_hash;
  }


};

