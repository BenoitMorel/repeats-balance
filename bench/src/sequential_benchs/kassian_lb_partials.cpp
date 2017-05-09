#include <iostream>
#include "common.h"


double get_unique_pattern_ratio(LikelihoodEngine &engine)
{
  double total_patterns = 0;
  double total_sites = 0;
  std::vector<Partition*> &partitions = engine.get_partitions();
  for (unsigned int p = 0; p < partitions.size(); ++p) {
    
    Partition *partition = partitions[p];
    total_sites += partition->get_sites_number();
    total_patterns += partition->get_unique_repeats_pattern_ratio() * double(partition->get_sites_number());
  }
  std::cout << "total sites " << total_sites << std::endl;
  return total_patterns / total_sites;
}

void treat_core(LikelihoodEngine &engine) 
{
  engine.update_partials();
}

void kassian_lb_partials(int argc, char *params[])
{
  if (argc != 8) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence partitions states use_repeats repeats_lookup_size iterations cores" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  const char *partition_file = params[i++];
  unsigned int states_number = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  unsigned int cores = atoi(params[i++]);

  
  MSA full_msa(seq, states_number); 
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(partition_file, initial_partitionning);
  std::vector<MSA *> msas;
  std::vector<WeightedMSA> weighted_msas;
  for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
    msas.push_back(new MSA(&full_msa, initial_partitionning[i], i));
    msas[i]->compress();
    weighted_msas.push_back(WeightedMSA(msas[i], 1.0));
  }
  LoadBalancer balancer;
  std::vector<CoreAssignment> assignments;
  balancer.kassian_load_balance(cores, weighted_msas, assignments);
  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  0, 
		  "avx"); 

  for (unsigned int core = 0; core < assignments.size(); ++core) {
    LikelihoodEngine engine(newick, msas, assignments[core], attribute, states_number, 4, repeats_lookup_size);
    engine.update_operations();
    engine.update_matrices();
    Timer timer;
    for (unsigned int i = 0; i < iterations; ++i) {
      treat_core(engine);
    }
    std::cout << "Pattern ratio: " << get_unique_pattern_ratio(engine) << std::endl;
    std::cout << timer.get_time() << "ms" << std::endl;
  }


}
