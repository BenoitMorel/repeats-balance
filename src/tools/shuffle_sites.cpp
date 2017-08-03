#include "tools.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include "../common/repeatsbalance.hpp"
  
void shuffle_sites(int argc, char *params[])
{
  if (argc != 3) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip partitions outputphylip" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  const char *outputphy_filename = params[i++];
 
  unsigned int states_number = 4;
  MSA msa(phy_filename, states_number);
  std::vector<PartitionIntervals> partitions_intervals;
  PartitionIntervals::parse(part_filename, partitions_intervals);
  std::vector<MSA *> submsas;
  std::vector< std::vector<unsigned int> > shuffled_indices(partitions_intervals.size());
  for (unsigned int i = 0; i < partitions_intervals.size(); ++i) 
  {
    submsas.push_back(new MSA(&msa, partitions_intervals[i], i));
    for (unsigned int j = 0; j < submsas[i]->get_length(); ++j) {
      shuffled_indices[i].push_back(j);
    }
    std::random_shuffle(shuffled_indices[i].begin(), shuffled_indices[i].end());
  }

  
  std::ofstream out(outputphy_filename);
  out << msa.get_tips_number() << " ";
  out << msa.get_length() << std::endl;
  for (unsigned int tip = 0; tip < msa.get_tips_number(); ++tip) {
    out << msa.get_pll_msa()->label[tip] << " "; 
    for (unsigned int p = 0; p < submsas.size(); ++p) {
      for (unsigned int site = 0; site < submsas[p]->get_length(); ++site) {
        unsigned int newsite = shuffled_indices[p][site];
        out << submsas[p]->get_pll_msa()->sequence[tip][newsite];
      }
    }
    out << std::endl;
  }
}



