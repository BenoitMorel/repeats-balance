#include "mpi_benchs.hpp"

void synchronized_ftt(int argc, char *params[])
{
  if (argc != 9) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "sequence partitions states use_repeats update_repeats repeats_lookup_size iterations randomized seed" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
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

  srand(seed); 
  MSA full_msa(seq, states_number);
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

  LikelihoodEngine engine(&tree, msas, assignments[rank_id], attribute, states_number, 4, repeats_lookup_size);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  Timer timer;
  for (unsigned int i = 0; i < iterations; ++i) {
    tree.randomize_pll_utree(&full_msa);
    engine.update_operations();
    engine.update_matrices();
    engine.update_partials(update_repeats);
  }
  if (!rank_id) {
    int local_time = timer.get_time();
    std::cout << rank_id << " " << local_time << "ms" << std::endl;
  }
  MPI_Finalize();
}
