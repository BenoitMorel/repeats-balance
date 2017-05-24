#include "mpi_benchs.hpp"

void synchronized_ftt(int argc, char *params[])
{
  if (argc != 13) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "sequence partitions states use_repeats update_repeats repeats_lookup_size iterations randomized seed "
      << "trees_number use_randomize_tree use_barrier use_update_operations"
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

  unsigned int trees_number = atoi(params[i++]);
  unsigned int use_randomize_tree = atoi(params[i++]);
  unsigned int use_barrier = atoi(params[i++]);
  unsigned int use_update_operations = atoi(params[i++]);
  
  if (!use_update_operations && ((trees_number > 1) || use_randomize_tree)) {
    std::cerr << "[ERROR] these parameter don't make any sense" << std::endl;
    return;
  }
  MPI_Init(NULL, NULL);

  int cores = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &cores);

  int rank_id = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);

  srand(seed); 

  if (!rank_id) {
    std::cout << "caption ";
    std::cout << "data: " << seq << ", ";
    std::cout << (randomized ? "lb_randomized" : "lb_kassian") << ", ";
    std::cout << (use_repeats ? "use_rep" : "no_use_rep") << ", ";
    std::cout << (update_repeats ? "update_rep" : "no_update_rep") << ", ";
    std::cout << "iterations: " << iterations << ", ";
    std::cout << "trees: " << trees_number << ", ";
    std::cout << (use_randomize_tree ? "shuffle_tree" : "no_shuffle_tree") << ", ";
    std::cout << (use_barrier ? "mpi_barriers" : "no_barriers") << ", ";
    std::cout << (use_update_operations ? "update_operations" : "no_update_operations") << ", ";
    std::cout << std::endl;
  }

  MSA full_msa(seq, states_number);
  std::vector<Tree *> trees;
  for (unsigned int i = 0; i < trees_number; ++i) {
      trees.push_back(new Tree(&full_msa));
  }
  Tree *current_tree = trees[0];
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(partition_file, initial_partitionning);
  std::vector<MSA *> msas;
  std::vector<WeightedMSA> weighted_msas;
  LoadBalancer balancer;
  for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
    msas.push_back(new MSA(&full_msa, initial_partitionning[i], i));
    msas[i]->compress();
  }

  Timer lb_timer;
  if (!randomized) {
    for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
      weighted_msas.push_back(WeightedMSA(msas[i], 1.0));
    }
  } else {
    std::cerr << rank_id << " before lb " << std::endl;
    balancer.compute_weighted_msa(msas, weighted_msas, PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX);
    std::cerr << rank_id << " after lb " << std::endl;
  }
  std::vector<CoreAssignment> assignments;
  balancer.kassian_load_balance(cores, weighted_msas, assignments);
  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  0, 
		  "avx");

  unsigned int lb_time = lb_timer.get_time();
  if (!rank_id) {
    std::cout << "Time spent in load balancing in rank " << rank_id << ": " << lb_time << "ms"  << std::endl;
  }
  LikelihoodEngine engine(current_tree, msas, assignments[rank_id], attribute, states_number, 4, repeats_lookup_size);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  MPI_Barrier(MPI_COMM_WORLD);
  Timer timer;
  for (unsigned int i = 0; i < iterations; ++i) {
    if (use_update_operations) {
      current_tree = trees[i%trees_number];
      if (use_randomize_tree) {
        current_tree->randomize_pll_utree(&full_msa);
      }
      engine.set_current_tree(current_tree);
      engine.update_operations();
    }
    engine.update_matrices();
    engine.update_partials(update_repeats);
    if (use_barrier) {
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  int local_time = timer.get_time();
  if (use_barrier) {
    if (!rank_id ) {
      std::cout << rank_id << " " << local_time << "ms" << std::endl;
    }
  } else {
    int max_time = 0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    std::cout << rank_id << " " << local_time << "ms" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank_id) {
      std::cout << "max_time: " << max_time << "ms" << std::endl;
    }
  }
  MPI_Finalize();

}
