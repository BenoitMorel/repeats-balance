#include "tools.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>
#include <math.h>

double get_distance(const std::vector<double> &v1, 
    const std::vector<double> &v2)
{
  double d = 0.0;
  for (unsigned int i = 0; i < v1.size(); ++i) {
    d += fabs(v1[i] * v1[i] - v2[i] * v2[i]);
  }
  return sqrt(d);
}

void add_vector(std::vector<double> &sum, 
    const std::vector<double> &v)
{
  for (unsigned int i = 0; i < v.size(); ++i) { 
    sum[i] += v[i];
  }
}

void multiply_vector(std::vector<double> &v, double factor) {
  for (unsigned int i = 0; i < v.size(); ++i) { 
    v[i] *= factor;
  }
}

void compute_weights(std::vector<double> &normalized_weights,
    const std::vector<Partition *> &partitions,
    Tree &tree)
{
  normalized_weights.resize(partitions.size());
  double sum = 0;
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_matrices(tree);
    partitions[i]->update_partials(tree);
    normalized_weights[i] = partitions[i] ->get_unique_repeats_pattern_ratio();
    sum += normalized_weights[i];
  }
  multiply_vector(normalized_weights, 1.0 / sum);
}

void display_weights_latex(std::vector<double> &weights) {
  for (unsigned int i = 0; i < weights.size(); ++i) {
    std::cout << "(" << i + 1 << "," << weights[i] << ")";
  }
  std::cout << std::endl;
}

void compute_subpartitions(std::vector<Partition *> &subpartitions,
    MSA &full_msa,
    const std::vector<PartitionIntervals> &partitions_intervals, 
    unsigned int attribute,
    unsigned int states_number)
{
  for (unsigned int i = 0; i < partitions_intervals.size(); ++i) {
    MSA submsa(&full_msa, partitions_intervals[i], i);
    submsa.compress();
    subpartitions.push_back(new Partition(&submsa, attribute, states_number, 4, 0));
  }
}

void analyse_partitions(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip partition states [newick|random trees number]" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  const char* tree_name = params[i++];

  // tree_name :
  // - integer n  -> look at n random trees
  // - newick file-> only look at this tree

  unsigned int tree_samples = 10;
  unsigned int trees_number = 1;
  const char *tree_filename = 0;
  if (atoi(tree_name)) {
    trees_number = atoi(tree_name);
  } else {
    tree_filename = tree_name;
  }

  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  srand(time(0));
  MSA msa(phy_filename, states_number);
  std::vector<PartitionIntervals> partitions_intervals;
  PartitionIntervals::parse(part_filename, partitions_intervals);
  unsigned int partitions_number = partitions_intervals.size();
  std::vector<Partition *> partitions;
  std::vector<double> temp_weights;
  std::vector<double> average_weights(partitions_number, 0.0);

  compute_subpartitions(partitions, msa, partitions_intervals, attribute, states_number);
  for (unsigned int i = 0; i < tree_samples; ++i) {
    Tree tree(&msa);
    tree.update_all_operations();
    compute_weights(temp_weights, partitions, tree);
    add_vector(average_weights, temp_weights);
  }
  multiply_vector(average_weights, 1.0 / tree_samples);
  display_weights_latex(average_weights);

  double max_distance = 0;
  std::vector<double> worst_weights;
  for (unsigned int i = 0; i < trees_number; ++i) {
    Tree tree(&msa);
    tree.update_all_operations();
    compute_weights(temp_weights, partitions, tree);
    double d = get_distance(average_weights, temp_weights);
    if (d > max_distance) {
      max_distance = d;
      std::cout << d << " ";
      tree.print();
      worst_weights = temp_weights;
    }
  }
  display_weights_latex(worst_weights);
}

