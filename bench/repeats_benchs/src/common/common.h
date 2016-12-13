#ifndef _LIBPLLBENCH_COMMON_H
#define _LIBPLLBENCH_COMMON_H

#include <time.h>
#include <vector>
#include "safepll.h"

struct PLLHelper {
  PLLHelper(const char *newick, const char *seq, unsigned int attribute);
  ~PLLHelper();
  void update_all_partials();
  void update_partials(int (*cbtrav)(pll_utree_t *), 
      unsigned int update_repeats, 
      unsigned int &traversal_size);

  void update_operations(int (*cbtrav)(pll_utree_t *), unsigned int &traversal_size);
  void update_partials(unsigned int update_repeats); 

  void update_partial(pll_operation_t *operation, 
      unsigned int iterations = 1, 
      bool update_repeats = true);
  void update_partials(std::vector<pll_operation_t> &operations, 
      unsigned int iterations = 1, 
      bool update_repeats = true);
  double get_likelihood();
  void save_svg(const char* file);
  void fill_children_number(std::vector<unsigned int> &children);
  // return true if path.size() == max_size
  bool build_path(unsigned int max_size, std::vector<pll_operation_t> &path);
  void build_random_path(double p, std::vector<pll_operation_t> &path);
  pll_utree_t *get_random_branch();

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


private:
  void fill_children_number_rec(pll_utree_t *root, std::vector<unsigned int> &children);
};

class Timer {
public:
  Timer();
  // ms
  long get_time();

private:
  timespec start;
};

#endif