#ifndef _TREE_
#define _TREE_

#include "safepll.h"
#include <vector>

class Tree {
public:
  Tree(const char *newick_file);
  
  virtual ~Tree();

  unsigned int get_tips_number() const {return tips_number;} 

  pll_utree_t *get_pll_tree() {return pll_utree;}
  
  const pll_utree_t *get_pll_tree() const {return pll_utree;}

  void update_operations(int (*traverse)(pll_utree_t *));
 
  const double *get_branch_lengths() const {return &branch_lengths[0];}

  const unsigned int *get_matrix_indices() const {return &matrix_indices[0];}

  unsigned int get_matrix_count() const {return matrix_indices.size();}

  const pll_operation_t *get_operations() const {return &operations[0];}
  
  unsigned int get_operations_number() const {return operations.size();}

private:
  static void set_missing_branch_length_recursive(pll_utree_t * tree, double length);
  static void set_missing_branch_length(pll_utree_t * tree, double length);

private:
  pll_utree_t * pll_utree;
  std::vector<double> branch_lengths;
  std::vector<unsigned int> matrix_indices;
  std::vector<pll_operation_t> operations;
  std::vector<pll_utree_t *> traverse_buffer;
  unsigned int tips_number;
  unsigned int innernodes_number;
  unsigned int nodes_number;
  unsigned int branches_number;
};

#endif
