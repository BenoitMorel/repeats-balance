
#include <iostream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"
#include "Helper.hpp"

void experience1() {
//  Helper::experiment1("../data/94/94.phy", "../data/94/94.part", 10, 30, "../results/experiment1/94_30cpus");
  Helper::experiment1("../data/94/94.phy", "../data/94/94.part", 10, 200, "../results/experiment1/94_200cpus");  
  Helper::experiment1("../data/59/59.phy", "../data/59/59.part", 10, 10, "../results/experiment1/59_10cpus");
  Helper::experiment1("../data/128/128.phy", "../data/128/128.part", 10, 10, "../results/experiment1/128_10cpus");
  Helper::experiment1("../data/404/404.phy", "../data/404/404.part", 10, 10, "../results/experiment1/404_10cpus");
}

void experience2() {
  Helper::experiment2("../data/59/59.phy", "../data/59/59.part", "../data/59/seed12345_top0.newick", 
                      10, 10, "../results/experiment2/59_10cpus_seed12345_top0");
  Helper::experiment2("../data/59/59.phy", "../data/59/59.part", "../data/59/seed12345_top79.newick", 
                      10, 10, "../results/experiment2/59_10cpus_seed12345_top79");
  Helper::experiment2("../data/128/128.phy", "../data/128/128.part", "../data/128/firsttree.newick", 
                      10, 10, "../results/experiment2/128_10cpus_first");
  Helper::experiment2("../data/128/128.phy", "../data/128/128.part", "../data/128/tree86894.newick", 
                      10, 10, "../results/experiment2/128_10cpus_tree86894");
  Helper::experiment2("../data/404/404.phy", "../data/404/404.part", "../data/404/seed42_firsttree.newick", 
                      10, 10, "../results/experiment2/404_10cpus_first");
  Helper::experiment2("../data/404/404.phy", "../data/404/404.part", "../data/404/seed42_tree_98240.newick", 
                      10, 10, "../results/experiment2/404_10cpus_98240");

}

void experience3() {
  Helper::experiment3("../data/128/split-partitions/12S_rRNA.phy", "../results/experiment3/exp3_12S_rRNA.tex");
  Helper::experiment3("../data/128/split-partitions/CREM.phy", "../results/experiment3/exp3_CREM.tex");
  Helper::experiment3("../data/128/split-partitions/BRCA1.phy", "../results/experiment3/exp3_BRCA1.tex");
}

int main()
{
  int seed = time(0);
  seed = 1477856748;
  std::cout << "seed : " << seed << std::endl;
  srand(seed);
  experience1();
  experience2(); 
  //experience3(); 
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
