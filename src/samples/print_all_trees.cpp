#include "samples.hpp"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "../common/repeatsbalance.hpp"


void print_rec(std::vector<unsigned int> &branch_seq, 
    unsigned int current_depth,
    MSA *msa) 
{
  if (current_depth == branch_seq.size()) {
    Tree tree(msa, branch_seq);
    tree.print();
  } else {
    for (unsigned int i = 0; i < 2 * current_depth - 3; ++i) {
      branch_seq[current_depth] = i;
      print_rec(branch_seq, current_depth + 1, msa);
    }
  }
}


void print_all_trees(int argc, char *params[])
{
  if (argc != 2) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip states" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);


  srand(time(0));
  MSA msa(phy_filename, states_number);
  std::vector<unsigned int> branch_seq(msa.get_tips_number());
  print_rec(branch_seq, 3, &msa); 
  

}

