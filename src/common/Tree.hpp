#ifndef _TREE_
#define _TREE_

#include "safepll.h"
#include <vector>

class MSA;

/**
 * Some of the code is copy-pasted from pll-modules
 */


class Tree {
public:
  // null newick_file for a random tree
  Tree(const MSA *msa, const char *newick_file = 0);
  
  // create a pll_utree_t creating a tree with the first 3 taxa in names, and
  // iteratively inserting the ith taxon from names at the branch branch_seq[i]
  // (3 first elements of branch_seq are ignored)
  Tree(const MSA *msa, const std::vector<unsigned int> &branch_seq);

  virtual ~Tree();

  unsigned int get_tips_number() const {return tips_number;} 

  pll_utree_t *get_pll_tree() {return pll_utree;}
  
  const pll_utree_t *get_pll_tree() const {return pll_utree;}

  pll_unode_t *get_pll_root() {return pll_utree->nodes[tips_number + innernodes_number - 1];}
  const pll_unode_t *get_pll_root() const {return pll_utree->nodes[tips_number + innernodes_number - 1];}

  void update_all_operations() {update_operations(traverser_full);}
  void update_operations(int (*traverse)(pll_unode_t *));
 
  const double *get_branch_lengths() const {return &branch_lengths[0];}

  const unsigned int *get_matrix_indices() const {return &matrix_indices[0];}

  unsigned int get_matrix_count() const {return matrix_indices.size();}

  const pll_operation_t *get_operations() const {return &operations[0];}
  
  unsigned int get_operations_number() const {return operations.size();}

  // update_operations must be called after that
  void randomize_pll_utree(const MSA *msa);

  static int traverser_full(pll_unode_t * node);

  char *get_newick() const;
  
  void print();
private:
  void init(const MSA *msa);
  
  pll_utree_t *create_from_vector(const std::vector<unsigned int> &branches_seq, 
      unsigned int taxa_count,
      const char * const* names);
  // create a random pll_utree_t of size taxa_count  
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
