#include <iostream>
#include "common.h"
#include <limits>



void naive_numerics(unsigned int iterations) {
  int ok_mult = 1;
  int ok_add = 1;
  for (unsigned int i = 0; i < iterations; ++i) {
    double a = ((double) rand() / (RAND_MAX)); 
    double b = ((double) rand() / (RAND_MAX));
    ok_mult &= (a*b*a*a*b*b) == (b*b*b*a*a*a);
    if (!ok_mult) {
      std::cout << (a*b*a*a*b*b) << " " << (b*b*b*a*a*a) << std::endl;
    }
    ok_add &= (a+b) == (b+a);
  }
  std::cout << ok_mult << " " << ok_add << std::endl;

}

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
void numerics(int argc, char *params[])
{
  if (argc != 6) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats iterations seed" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  unsigned int seed = atoi(params[i++]);

  srand(seed);
  std::cout.precision(std::numeric_limits< double >::digits10);
  std::cerr.precision(std::numeric_limits< double >::digits10);
  const char *arch="avx";
  unsigned int additional_attr=0;

  unsigned int attribute = PLLHelper::compute_attribute(use_repeats, additional_attr, arch);
  if (INVALID_ATTRIBUTE == attribute) {
    return;
  }
  PLLHelper helper(newick, seq, attribute, states); 
  
 
  for (unsigned int i = 0; i < iterations; ++i) {
    helper.generate_random_model();
    helper.update_all_partials();
    std::cerr << "l " << helper.get_likelihood() << std::endl;
    double d[2];
    helper.get_derivative(d , &d[1]); 
    std::cerr << "d " << d[0] << " " << d[1]  << std::endl;
  }


}



