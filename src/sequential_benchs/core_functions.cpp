#include <iostream>
#include "common.h"

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
void core_functions(int argc, char *params[])
{
  if (argc != 10) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats  additional_attr repeats_lookup_size "
      << "iterations_likelihood iterations_sumtable iterations_derivatives  arch"
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
  unsigned int iterations_likelihood = atoi(params[i++]);
  unsigned int iterations_sumtable = atoi(params[i++]);
  unsigned int iterations_derivatives = atoi(params[i++]);
  const char *arch = params[i++];

  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  additional_attr, 
		  arch); 
  LikelihoodEngine engine(newick, seq, 0, attribute, states, 4, repeats_lookup_size);
  
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  double ll = engine.compute_likelihood();
  double d_f = 0;
  double dd_f = 0;
  engine.update_sumtable();
  engine.compute_derivatives(&d_f, &dd_f);
 

  Timer timer;
  for (i = 0; i < iterations_likelihood; ++i) {
    ll = engine.compute_likelihood();
  } 
  for (i = 0; i < iterations_sumtable; ++i) {
    engine.update_sumtable();
  } 
  for (i = 0; i < iterations_derivatives; ++i) {
    engine.compute_derivatives(&d_f, &dd_f);
  } 
  std::cerr << "ll = " << ll << " derivatives " <<  d_f << " " <<  dd_f << std::endl; 
  std::cout << " " << timer.get_time() << "ms" << std::endl; 
}



