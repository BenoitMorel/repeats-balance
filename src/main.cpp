
#include <iostream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"


void compute_average_SRcount(const std::string &sequences_file_name) {
  
  InputSequences sequences;
  const unsigned int iterations = 100;
  parse_sequences(sequences_file_name.c_str(), sequences);  
  SRLOG("conpute_average_SRcount with " << sequences.number() << " seq of size " << sequences.width());
  Tree tree;
  std::vector<int> SRCount(sequences.width());
  for (unsigned int i = 0; i < iterations; ++i) {
    SRLOG("iteration" << i);
    tree.set_random(sequences.number(), time(0) + i);
    tree.update_SRcount(sequences, SRCount); 
  }
  //  for (unsigned int site = 0; site < SRCount.size(); ++site) {
  //  std::cout << (double)SRCount[site] / ((double)iterations * (sequences._seq_number - 1)) << " ";
  //}  
  std::cout << std::endl;
}

int main()
{
  //test1();
  //test2();
  //test_sequences_parser("../data/minimal-6/minimal-6.phy");
  //test_print_random_trees();
  //compute_average_SRcount("../data/minimal-6/minimal-6.phy");
  //compute_average_SRcount("../data/simple_seq/simple4-4.phy");
  compute_average_SRcount("../data/128/128.phy");
  return 0;
}
