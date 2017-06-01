#include <iostream>
#include "common.h"

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
void full_traversal(int argc, char *params[])
{
  if (argc != 8) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats additional_attr repeats_lookup_size iterations arch" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int additional_attr = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  const char *arch = params[i++];


  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  additional_attr, 
		  arch); 
  LikelihoodEngine engine(newick, seq, 0, attribute, states, 4, repeats_lookup_size);
  
  engine.update_operations();
  engine.update_matrices();
  Timer timer;
  for (i = 0; i < iterations; ++i) {
    engine.update_partials();
  } 
  std::cout << " " << timer.get_time() << "ms" << std::endl; 
  std::cerr << engine.compute_likelihood() << std::endl; 
}


