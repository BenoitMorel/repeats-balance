#include "tools.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
  
void generate_partitions(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip partitions_number model output_part" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  unsigned int partitions_number = atoi(params[i++]);
  const char *model = params[i++];
  const char *part_filename = params[i++];

  std::ifstream phy_file(phy_filename);
  std::ofstream part_file(part_filename);
  if (!phy_file.is_open()) {
    std::cout << "Cannot open : " << phy_filename << std::endl; 
    return;
  }
  if (!part_file.is_open()) {
    std::cout << "Cannot open : " << part_filename << std::endl; 
    return;
  }
  unsigned int sites = 0;
  phy_file  >> sites;
  phy_file  >> sites;
  unsigned int length = sites / partitions_number + 1;
  for (unsigned int s = 0; s < sites; s += length) {
    part_file << model << ", Partition" << s/length << " = " << s + 1 << "-" << std::min(sites, s + length) << std::endl; 
  }
}



