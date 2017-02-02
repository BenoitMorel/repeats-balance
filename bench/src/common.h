#ifndef _LIBPLLBENCH_COMMON_H
#define _LIBPLLBENCH_COMMON_H

#include <time.h>
#include <vector>
#include "safepll.h"

#ifndef PLL_ATTRIB_SITES_REPEATS // old pll version
  #define PLL_ATTRIB_SITES_REPEATS 0
#else
  #define HAS_REPEATS
#endif

#define INVALID_ATTRIBUTE ((unsigned int)1)

struct PLLHelper {
  static unsigned int compute_attribute(bool use_repeats, 
                                        unsigned int additional_attr, 
                                        const char *arch);
  PLLHelper(const char *newick, const char *seq, unsigned int attribute, unsigned int states = 4); 
  ~PLLHelper();
  void set_srlookup_size(unsigned int size);
  
  void update_all_partials(bool update_repeats = true);
  void update_partials(int (*cbtrav)(pll_utree_t *), 
      unsigned int update_repeats, 
      unsigned int &traversal_size);
  void update_all_repeats();
  void update_partials(unsigned int update_repeats); 
  void update_partial(pll_operation_t &operation, 
      unsigned int iterations = 1, 
      bool update_repeats = true);

  
  void update_operations(int (*cbtrav)(pll_utree_t *), unsigned int &traversal_size);
  double get_likelihood();
  void update_sumtables(); 
  void get_derivative(double *d_f, double *dd_f); 

  // call update_all_partials at least once before
  void init_tree_stats();
  // prints : op_left_children_number, op_right_children_number, op_depth
  void print_op_stats(pll_operation_t &op) const;

  // class members
  pll_utree_t * tree;
  pll_partition_t *partition;
  double * branch_lengths;
  unsigned int * matrix_indices;
  unsigned int params_indices[4];
  unsigned int tip_nodes_count;
  unsigned int inner_nodes_count;
  unsigned int nodes_count;
  unsigned int branch_count;
  unsigned int ops_count;
  unsigned int matrix_count;
  pll_utree_t ** travbuffer;
  pll_operation_t * operations;
  double * sumtable;
  // these values are indexed by the clv index
  std::vector<unsigned int> srclasses_number;
  std::vector<unsigned int> children_number;
  std::vector<unsigned int> depths;



private:
  // compute children_number and depths
  void init_tree_stats_rec(pll_utree_t *root, unsigned int depth);
  void update_partials(std::vector<pll_operation_t> &operations, 
      unsigned int iterations = 1, 
      bool update_repeats = true);
};

class Timer {
public:
  Timer();
  // ms
  long get_time();

private:
  timespec start;
};


void full_traversal(int argc, char *params[]);
void update_repeats(int argc, char *params[]);
void derivatives(int argc, char *params[]);
void sumtables(int argc, char *params[]);
void partitioned_full_traversal(int argc, char *params[]);
void pernode(int argc, char *params[]);

#endif
