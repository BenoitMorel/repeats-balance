#include "Tree.hpp"
#include <iostream>


Tree::Tree(const char *newick_file)
{
  unsigned int tips_number = 0;
  pll_utree = pll_utree_parse_newick(newick_file, &tips_number);
  if (!pll_utree) {
    std::cerr << "Error: Tree::tree null pll_utree" << std::endl;
  }
  set_missing_branch_length(pll_utree, 0.000001);

  innernodes_number = tips_number - 2;
  nodes_number = innernodes_number + tips_number;
  branches_number = nodes_number - 1;
  
  branch_lengths.resize(branches_number);
  matrix_indices.resize(branches_number);

}



Tree::~Tree()
{

}

void Tree::update_operations(int (*traverse)(pll_utree_t *))
{
  unsigned int traverse_size = 0;
  traverse_buffer.resize(nodes_number);
  pll_utree_traverse(pll_utree,
    traverse,
    &traverse_buffer[0],
    &traverse_size);
  traverse_buffer.resize(traverse_size);

  operations.resize(branches_number);
  matrix_indices.resize(branches_number);

  unsigned int matrices_number = 0;
  unsigned int operations_number = 0;
  pll_utree_create_operations(&traverse_buffer[0],
      traverse_size,
      &branch_lengths[0],
      &matrix_indices[0],
      &operations[0],
      &matrices_number,
      &operations_number);

  matrix_indices.resize(matrices_number);
  operations.resize(operations_number);
}

void Tree::set_missing_branch_length_recursive(pll_utree_t * tree, double length)
{
  if (tree) {
    /* set branch length to default if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next) {   
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive(tree->next->back, length);
        set_missing_branch_length_recursive(tree->next->next->back, length);
    }   
  }
}

void Tree::set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
}

