
#include <iostream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"
#include "Helper.hpp"
int main()
{
  srand(time(0));
  //compute_average_SRcount("../data/simple_seq/simple4-4.phy");
  //compute_partitions_sr_count("../data/128/128.phy", "../data/128/128.part", 50);
  Helper::print_stats("../data/128/128.phy", "../data/128/128.part", 1, 5, "stats128.txt");
  Helper::count_sr("../data/128/128.phy", "../data/128/single_partition.part","../data/128/RAxML_parsimonyTree.128");
  Helper::count_sr("../data/128/128.phy", "../data/128/subdivisions.part","../data/128/RAxML_parsimonyTree.128");
  Helper::count_sr("../data/libplldata/unrooted.phy", "../data/libplldata/single.part","../data/libplldata/rooted.newick");
  Helper::count_sr("../data/libplldata/plop.phy", "../data/libplldata/single.part","../data/libplldata/rooted.newick");
  
  return 0;
}


void compute_average_SRcount(const std::string &sequences_file_name) {
  
  InputSequences sequences;
  const unsigned int iterations = 1;
  parse_sequences(sequences_file_name.c_str(), sequences);  
  SRLOG("conpute_average_SRcount with " << sequences.number() << " seq of size " << sequences.width());
  Tree tree;
  std::vector<double> SRCount(sequences.width());
  for (unsigned int i = 0; i < iterations; ++i) {
    SRLOG("iteration" << i);
    tree.set_random(sequences.number());
    tree.update_SRcount(sequences, 0, SRCount);
	std::cout << tree << std::endl; 
  }
  for (unsigned int site = 0; site < SRCount.size(); ++site) {
    std::cout << (double)SRCount[site] / ((double)iterations * (sequences.number() - 1)) << " ";
  }  
  std::cout << std::endl;
}

void compute_partitions_sr_count(const std::string &sequences_file_name, const std::string &partitions_file_name, unsigned int iterations) {
  InputSequences sequences;
  parse_sequences(sequences_file_name.c_str(), sequences);  
  InputPartitions inputpartitions;
  parse_partitions(partitions_file_name.c_str(), inputpartitions);  
  Partitions partitions;
  inputpartitions.generate_partitions(partitions, &sequences);

  std::cout << "compute_partitions_sr_count" << std::endl;
  std::cout << "-- sequences number :     " << sequences.number() << std::endl;
  std::cout << "-- total sequences size : " << sequences.width() << std::endl;
  std::cout << "-- partitions number :    " << partitions.size() << std::endl;
  Tree tree;
  for (unsigned int i = 0; i < iterations; ++i) {
    SRLOG("iteration" << i);
    tree.set_random(sequences.number());
    for (unsigned int j = 0; j < partitions.size(); ++j) {
      tree.update_SRcount(sequences, partitions[j].start(), partitions[j].site_costs()); 
    }
  }
  // max site repeats count : iterations * number of inner nodes
  for (unsigned int j = 0; j < partitions.size(); ++j) {
    partitions[j].normalize_costs(double(iterations * (sequences.number() - 1)));
  }
  for (unsigned int j = 0; j < partitions.size(); ++j) {
    std::cout << partitions[j] << std::endl;
  }
}
