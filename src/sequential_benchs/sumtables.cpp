#include <iostream>
#include "common.h"

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
void sumtables(int argc, char *params[])
{
  if (argc != 9) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats update_repeats additional_attr repeats_lookup_size iterations arch" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  i++;
  unsigned int additional_attr = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  const char *arch = params[i++];


  unsigned int attribute = PLLHelper::compute_attribute(use_repeats, additional_attr, arch);
  if (INVALID_ATTRIBUTE == attribute) {
    return;
  }
  PLLHelper helper(newick, seq, attribute, states); 
  std::cout << "srlookup" << std::endl;
  helper.set_srlookup_size(repeats_lookup_size);
  std::cout << "update all" << std::endl;
  helper.update_all_partials();
  std::cout << "done" << std::endl;
  Timer t;
  double d_f, dd_f;
  for (i = 0; i < iterations; ++i) {
    helper.update_sumtables();  
  }
  helper.get_derivative(&d_f, &dd_f);
  std::cerr << "derivatives " <<  d_f << " " <<  dd_f << std::endl; 
  std::cout << " " << t.get_time() << "ms" << std::endl; 
}


