#ifndef _TREE_
#define _TREE_

#include "safepll.h"
#include <vector>

class MSA;

class Tree {
public:
  // null newick_file for a random tree
  Tree(const MSA *msa, const char *newick_file = 0);

  virtual ~Tree();

  unsigned int get_tips_number() const {return tips_number;} 

  pll_utree_t *get_pll_tree() {return pll_utree;}
  
  const pll_utree_t *get_pll_tree() const {return pll_utree;}

  const pll_unode_t *get_pll_root() const {return pll_utree->nodes[tips_number + innernodes_number - 1];}

  void update_operations(int (*traverse)(pll_unode_t *));
 
  const double *get_branch_lengths() const {return &branch_lengths[0];}

  const unsigned int *get_matrix_indices() const {return &matrix_indices[0];}

  unsigned int get_matrix_count() const {return matrix_indices.size();}

  const pll_operation_t *get_operations() const {return &operations[0];}
  
  unsigned int get_operations_number() const {return operations.size();}

  // update_operations must be called after that
  void randomize_pll_utree(const MSA *msa);

  static int traverser_full(pll_unode_t * node);

private:
  void init(const MSA *msa);
  pll_utree_t *create_random(unsigned int taxa_count, const char * const* names);
  void map_to_labels(const MSA* msa);
private:
  pll_utree_t * pll_utree;
  std::vector<double> branch_lengths;
  std::vector<unsigned int> matrix_indices;
  std::vector<pll_operation_t> operations;
  std::vector<pll_unode_t *> traverse_buffer;
  unsigned int tips_number;
  unsigned int innernodes_number;
  unsigned int nodes_number;
  unsigned int branches_number;
};

#endif
