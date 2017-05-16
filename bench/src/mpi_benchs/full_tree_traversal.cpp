#include "mpi_benchs.hpp"

double print_stats(LikelihoodEngine &engine)
{
  double total_patterns = 0;
  double total_sites = 0;
  std::vector<Partition*> &partitions = engine.get_partitions();
  for (unsigned int p = 0; p < partitions.size(); ++p) {
    
    Partition *partition = partitions[p];
    total_sites += partition->get_sites_number();
    total_patterns += partition->get_unique_repeats_pattern_ratio() * double(partition->get_sites_number());
  }
  double res = total_patterns / total_sites;
  std::cout << "  total sites " << total_sites << std::endl;
  std::cout << "  ratio " << res << std::endl;
  std::cout << "  Sites * ratio = " << total_patterns << std::endl;
  return res;
}

void full_tree_traversal(int argc, char *params[])
{
  if (argc != 10) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence partitions states use_repeats update_repeats repeats_lookup_size iterations randomized seed" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  const char *partition_file = params[i++];
  unsigned int states_number = atoi(params[i++]);
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int update_repeats = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  unsigned int randomized = atoi(params[i++]);
  unsigned int seed = atoi(params[i++]);


  MPI_Init(NULL, NULL);

  int cores = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &cores);

  int rank_id = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);

  if (!rank_id) {
    std::cout << "caption ";
    std::cout << cores << "cores ";
    std::cout << (use_repeats ? "repeats " : "tipinner ");
    std::cout << (update_repeats ? "update " : "noupdate ")  << std::endl;
    std::cout << (randomized ? "randomized " : "kassian ")  << std::endl;
  }
  
    
  srand(seed); 
  MSA full_msa(seq, states_number); 
  //Tree tree(newick);
  Tree tree(&full_msa);
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(partition_file, initial_partitionning);
  std::vector<MSA *> msas;
  std::vector<WeightedMSA> weighted_msas;
  LoadBalancer balancer;
  for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
    msas.push_back(new MSA(&full_msa, initial_partitionning[i], i));
    msas[i]->compress();
  }
  if (!randomized) {
    for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
      weighted_msas.push_back(WeightedMSA(msas[i], 1.0));
    }
  } else {
    balancer.compute_weighted_msa(msas, weighted_msas, PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX);
  }
  std::vector<CoreAssignment> assignments;
  balancer.kassian_load_balance(cores, weighted_msas, assignments);
  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  0, 
		  "avx");

  if (!rank_id) {
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      std::cout << assignments[i] << std::endl;
    }
  }
  LikelihoodEngine engine(&tree, msas, assignments[rank_id], attribute, states_number, 4, repeats_lookup_size);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  Timer timer;
  for (unsigned int i = 0; i < iterations; ++i) {
    engine.update_partials(update_repeats);
  }
  int local_time = timer.get_time();
  print_stats(engine);
  std::cout << rank_id << " " << local_time << "ms" << std::endl;
  MPI_Finalize();
}
