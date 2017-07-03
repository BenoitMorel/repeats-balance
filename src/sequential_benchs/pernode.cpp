#include <iostream>
#include "common.h"

void print_op_stats(const Partition &partition, const Tree &tree, unsigned int op)
{
  std::cout << op << " ";
  pll_operation_t operation = tree.get_operations()[op];
  std::cout << pll_get_sites_number(partition.get_partition(), operation.parent_clv_index) << " ";
  std::cout << pll_get_sites_number(partition.get_partition(), operation.child1_clv_index) << " ";
  std::cout << pll_get_sites_number(partition.get_partition(), operation.child2_clv_index);
  std::cout << std::endl;
}

void update_partial(Partition &partition, 
    Tree &tree, 
    unsigned int op,
    unsigned int update_repeats,
    unsigned int iterations)
{
  std::vector<pll_operation_t> operations(iterations, tree.get_operations()[op]);
  partition.update_partials_from_operations(&operations[0], iterations, update_repeats); 
}

/*
 *  For each node, print : left and right children number, and clv_computation time
 */
void pernode(int argc, char *params[])
{
  if (argc != 10 && argc != 7) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states arch node_index repeats_lookup_size additional_attr  use_repeats update_repeats iterations" 
      << std::endl;
    std::cerr 
      << "newick sequence states arch node_index repeats_lookup_size additional_attr" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  const char *arch = params[i++];
  int node_index =  atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int additional_attr = atoi(params[i++]);

  MSA msa(seq, states);
  msa.compress();
  Tree tree(&msa, newick);
  tree.update_all_operations();
  
    // just print nodes stats
  if (argc == 7) {
    unsigned int attribute =  Partition::compute_attribute(1, 0, arch);
    Partition partition(&msa, attribute, states, 4,  repeats_lookup_size); 
    partition.update_matrices(tree);
    partition.update_partials(tree);  
    if (node_index == -1) {
      for (unsigned int op = 0; op < tree.get_operations_number(); ++op) {
        print_op_stats(partition, tree, op);
      }
    } else {
      print_op_stats(partition, tree, node_index);
    }
    return;
  }
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int update_repeats = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  unsigned int attribute = Partition::compute_attribute(use_repeats, additional_attr, arch);
  Partition partition(&msa, attribute, states, 4,  repeats_lookup_size); 
  partition.update_partials(tree);  
  
  if (node_index == -1) {
    for (unsigned int op = 0; op < tree.get_operations_number(); ++op) {
      Timer t;
      update_partial(partition, tree, op, update_repeats, iterations);
      std::cout << t.get_time() << "ms" << std::endl; 
    }
  } else {
    Timer t;
    update_partial(partition, tree, node_index, update_repeats, iterations);
    std::cout << t.get_time() << "ms" << std::endl; 
  }
}



