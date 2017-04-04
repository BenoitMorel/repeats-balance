#include <iostream>
#include "common.h"
#include <fstream>

void export_operations(pll_operation_t * operations,
    unsigned int operations_number,
    const char *output_file)
{
  std::cout << "Exporting operations into " << output_file << std::endl;
  std::ofstream os(output_file);
  os << operations_number << std::endl;
  for (unsigned int i = 0; i < operations_number; ++i) {
    os << operations[i].parent_clv_index << " ";
    os << operations[i].child1_clv_index << " ";
    os << operations[i].child2_clv_index;
    os << std::endl;
  }
}

void export_repeats_partition(pll_partition_t *partition,
    const char* output_file)
{
  std::cout << "Exporting partition into " << output_file << std::endl;
  std::ofstream os(output_file);
  os << partition->sites << std::endl;
  unsigned int nodes = partition->tips + partition->clv_buffers;
  for (unsigned int n = 0; n < nodes; ++n) {
    unsigned int *site_id = partition->repeats->pernode_site_id[n];
    for (unsigned int s = 0; s < partition->sites; ++s) {
      os << (partition->repeats->pernode_max_id[n] ? site_id[s] : s) << " ";
    }
    os << std::endl;
  }
}

/*
 *  export the repeats structure for the hypergraph project
 */
void export_repeats(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick split_sequences part_numbers output_dir" 
      << std::endl;
    std::cerr << "split_sequences is the directories containing PARTITION_x.phy files"<< std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seqdir = params[i++];
  unsigned int part_number = atoi(params[i++]);
  const char *output_dir = params[i++];

  char seq[1000];
  char output_file[1000];
  for(unsigned int i = 0; i < part_number; ++i) {
    sprintf(seq, "%sPARTITION_%d.phy", seqdir, i);
    PLLHelper helper(newick, seq, PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX); 
    helper.update_all_partials();
    if (i == 0) {
      sprintf(output_file, "%s/operations.txt", output_dir);
      export_operations(helper.operations, helper.ops_count, output_file);
    }
    sprintf(output_file, "%s/repeats_p%d.txt", output_dir, i);
    export_repeats_partition(helper.partition, output_file); 
  }
}





