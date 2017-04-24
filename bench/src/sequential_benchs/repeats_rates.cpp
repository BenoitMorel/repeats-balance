#include <iostream>
#include "common.h"

/*
 *  For each node, print : left and right children number, and clv_computation time
 */
void repeats_rates(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick split_sequences part_numbers states" 
      << std::endl;
    std::cerr << "split_sequences is the directories containing PARTITION_x.phy files"<< std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seqdir = params[i++];
  unsigned int part_number = atoi(params[i++]);

  char seq[1000];
  for(unsigned int i = 0; i < part_number; ++i) {
    sprintf(seq, "%sPARTITION_%d.phy", seqdir, i);
    PLLHelper helper(newick, seq, PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX); 
    helper.update_all_partials();
    std::cout << "p" << i << " " << helper.get_repeats_rates() * 100 << std::endl;
  }
  std::cout << std::endl;
}




