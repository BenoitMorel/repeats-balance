#include "tools.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>

void analyse_partitions(int argc, char *params[])
{
  if (argc != 3 && argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip partition states (tree)" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  const char* newick_filename = (argc == 4) ? params[i++] : 0;

  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  srand(time(0));
  MSA msa(phy_filename, states_number);
  std::vector<PartitionIntervals> partitions_intervals;
  PartitionIntervals::parse(part_filename, partitions_intervals);
  Tree tree(&msa, newick_filename);
  tree.update_operations(Tree::traverser_full);
  std::vector<double> pattern_compressions;
  std::vector<unsigned int> compressed_sites;
  std::vector<double> partition_weights;
  double weights_sum = 0;
  for (unsigned int i = 0; i < partitions_intervals.size(); ++i) {
    MSA submsa(&msa, partitions_intervals[i], i);
    unsigned int initial_sites = submsa.get_pll_msa()->length;
    submsa.compress();
    Partition partition(&submsa, attribute, states_number, 4, 0);
    partition.update_matrices(tree);
    partition.update_partials(tree);
    pattern_compressions.push_back( double(partition.get_sites_number())/double(initial_sites));
    compressed_sites.push_back(partition.get_sites_number());
    partition_weights.push_back(partition.get_unique_repeats_pattern_ratio());
    weights_sum += partition_weights[partition_weights.size() - 1];
  }
  for (unsigned int i = 0; i < partitions_intervals.size(); ++i) {
    std::cout << "Partition " << i << ":" << std::endl;
    std::cout << "  pattern compression " << pattern_compressions[i] << std::endl;
    std::cout << "  compressed sites: " << compressed_sites[i] << std::endl;
    std::cout << "  partition_weight: " << partition_weights[i] << std::endl;
    std::cout << "  normalized weights: " << partition_weights[i] / weights_sum << std::endl;

  }
}

