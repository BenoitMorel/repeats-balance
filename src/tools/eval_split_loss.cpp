#include "tools.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>
#include <math.h>


double get_repeats_rate(const MSA &msa, Tree &tree, unsigned int attribute) 
{
  Partition partition(&msa, attribute, msa.get_states_number(),
    4, 0); 
  partition.update_matrices(tree);
  partition.update_partials(tree);
  return partition.get_unique_repeats_pattern_ratio();
}

void split_msa(const MSA &msa, unsigned int splits_number,
	std::vector<MSA *> &output_msas)
{
  unsigned int size = msa.get_length() / splits_number + 1;
  for (unsigned int i = 0; i < msa.get_length(); i += size) {
    PartitionIntervals interval(0);
    unsigned int s = std::min(size, msa.get_length() - i);
    interval.add_interval(i, s);
    std::cout << s << " " << msa.get_length() << std::endl;
    output_msas.push_back(new MSA(&msa, interval, 0));
  }
}



void eval_split_loss(int argc, char *params[])
{
  if (argc != 3) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip states splits_number" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  unsigned int splits_number = atoi(params[i++]);
  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  srand(time(0));
  MSA msa(phy_filename, states_number);
  msa.compress();
  Tree tree(&msa);
  tree.update_all_operations();
  std::vector<MSA *> split_msas;
  split_msa(msa, splits_number, split_msas);
  std::cout << "Initial repeats rate: " << 
    get_repeats_rate(msa, tree, attribute)<< std::endl;;
  double total_rate = 0.0;
  for (unsigned int i = 0; i < split_msas.size(); ++i) {
    double rate = get_repeats_rate(*split_msas[i], tree, attribute);
    total_rate += rate / splits_number;
    std::cout << rate << std::endl;
  }
  std::cout << "Final repeats rate: " << total_rate << std::endl;

}

