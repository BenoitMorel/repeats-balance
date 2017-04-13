
#include "Tree.hpp"
#include "Partition.hpp"
#include "LikelihoodEngine.hpp"
#include <iostream>

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
int main(int argc, char *params[])
{
  if (argc != 10) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats update_repeats additional_attr repeats_lookup_size iterations arch" 
      << std::endl;
    return 1;
  }
  unsigned int i = 1;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int update_repeats = atoi(params[i++]);
  unsigned int additional_attr = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  const char *arch = params[i++];

  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  additional_attr, 
		  arch); 
  Tree tree(newick);
  Partition partition(seq, tree, attribute, states);
  LikelihoodEngine engine(tree, partition);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  std::cout << "Likelihood : " << engine.compute_likelihood() << std::endl;


  return 0;
}


