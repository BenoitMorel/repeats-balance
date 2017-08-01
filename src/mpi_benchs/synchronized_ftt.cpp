#include "mpi_benchs.hpp"

void synchronized_ftt(int argc, char *params[])
{
  if (argc != 12) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "sequence partitions states use_repeats update_repeats repeats_lookup_size iterations randomized seed "
      << "sample_tree trees_to_traverse use_barrier"
      << std::endl;
    std::cerr << "sample_tree can be a newick file or the number of random trees to generate" << std::endl;
    std::cerr << "trees_to_traverse can be a newick file or \"random\'" << std::endl;
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
  const char*  sample_tree_param = params[i++];
  const char*  newick_tree = params[i++];
  unsigned int use_barrier = atoi(params[i++]);
  
  MPI_Init(NULL, NULL);

  int cores = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &cores);

  int rank_id = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);

  srand(seed); 

  if (!rank_id) {
    for(i = 0; i < argc; i++) {
      printf("%s ", params[i]);
    }
    printf("\n");
    std::cout << "caption ";
    std::cout << "data: " << seq << ", ";
    std::cout << "part: " << partition_file << ", ";
    std::cout << (randomized ? "lb_randomized" : "lb_kassian") << ", ";
    std::cout << (use_repeats ? "use_rep" : "no_use_rep") << ", ";
    std::cout << (update_repeats ? "update_rep" : "no_update_rep") << ", ";
    std::cout << "iterations: " << iterations << ", ";
    std::cout << "tree: " << newick_tree << ", ";
    std::cout << (use_barrier ? "mpi_barriers" : "no_barriers") << ", ";
    std::cout << std::endl;
  }


  // ********************
  //  INIT
  // ********************
  MSA full_msa(seq, states_number);
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(partition_file, initial_partitionning);
  std::vector<MSA *> msas;
  std::vector<WeightedMSA> weighted_msas;
  LoadBalancer balancer;
  for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
    msas.push_back(new MSA(&full_msa, initial_partitionning[i], i));
    msas[i]->compress();
  }

  std::vector<Tree *> trees_sample;
  if (atoi(sample_tree_param)) {
    unsigned int sample_size = atoi(sample_tree_param);
    for (unsigned int i = 0; i < sample_size; ++i) {
      trees_sample.push_back(new Tree(&full_msa));
    }
  } else {
    trees_sample.push_back(new Tree(&full_msa, sample_tree_param));
  }
  for (unsigned int i = 0; i < trees_sample.size(); ++i) {
    trees_sample[i]->update_all_operations();
  }
  
  // ********************
  //  LOAD BALANCE
  // ********************
  Timer lb_timer;
  if (!randomized) {
    for (unsigned int i = 0; i < initial_partitionning.size(); ++i) {
      weighted_msas.push_back(WeightedMSA(msas[i], 1.0));
    }
  } else {
    balancer.compute_weighted_msa(msas, weighted_msas, PLL_ATTRIB_SITE_REPEATS | PLL_ATTRIB_ARCH_AVX, trees_sample);
  }
  std::vector<CoreAssignment> assignments;
  balancer.kassian_load_balance(cores, weighted_msas, assignments);
  if (rank_id == 0) {
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      std::cout << assignments[i] << std::endl;
    }
  }
  unsigned int lb_time = lb_timer.get_time();
  if (!rank_id) {
    std::cout << "Time spent in load balancing in rank " << rank_id << ": " << lb_time << "ms"  << std::endl;
  }
  // ********************
  //  BENCH
  // ********************
  newick_tree = memcmp("random", newick_tree, strlen("random")) ? newick_tree : 0;
  bool use_random_trees = !newick_tree;
  Tree current_tree(&full_msa, newick_tree);
  unsigned int attribute = Partition::compute_attribute(use_repeats, 
		  0, 
		  "avx");
  //attribute = attribute | PLL_ATTRIB_NOBCLV;
  LikelihoodEngine engine(&current_tree, msas, assignments[rank_id], attribute, states_number, 4, repeats_lookup_size);
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  MPI_Barrier(MPI_COMM_WORLD);
  Timer timer;
  for (unsigned int i = 0; i < iterations; ++i) {
    if (!use_random_trees) {
      current_tree.randomize_pll_utree(&full_msa);
    }
    //engine.set_current_tree(current_tree);
    engine.update_operations();
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
