#include <iostream>
#include "common.h"
#include <fstream>

void export_operations(const pll_operation_t * operations,
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
      << "newick phy part_file output_dir" 
      << std::endl;
    std::cerr << "split_sequences is the directories containing PARTITION_x.phy files"<< std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick_file = params[i++];
  const char *seq_file = params[i++];
  const char *part_file = params[i++];
  const char *output_dir = params[i++];
  unsigned int states = 4;
  unsigned int repeats_lookup_size = 2000000;
  LikelihoodEngine engine(newick_file, seq_file, part_file,
      PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX, states, 4, repeats_lookup_size);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();

  char output_file[1000];
  sprintf(output_file, "%s/operations.txt", output_dir);
  export_operations(engine.get_tree().get_operations(), engine.get_tree().get_operations_number(), output_file);

  std::vector<Partition*> partitions = engine.get_partitions();
  for(unsigned int i = 0; i < partitions.size(); ++i) {
    sprintf(output_file, "%s/repeats_p%d.txt", output_dir, i);
    export_repeats_partition(partitions[i]->get_partition(), output_file); 
  }
}





