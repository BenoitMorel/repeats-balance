#ifndef _HELPER_H_
#define _HELPER_H_

#include <string>
#include <fstream>
#include <vector>
#include "repeatsbalance.h"
#include "Partition.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"

class Helper {
public:
  static void print_count_sr(const std::string &sequences_file,
                          const std::string &partitions_file,
                          const std::string &tree_file
                          ) {

    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences);
    InputPartitions inputPartitions;
    parse_partitions(partitions_file.c_str(), inputPartitions);
    Partitions partitions;
    inputPartitions.generate_partitions(partitions, &sequences);
    Tree tree;
    parse_tree(tree_file.c_str(), sequences, tree);
    std::cout << "Average site repeat rate for each partition : " << std::endl;
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      tree.update_SRcount(partitions[i]);
      partitions[i].normalize_costs(sequences.number() - 1);
      std::cout << inputPartitions.name(i) << ":" << partitions[i].average_sr_rate() << std::endl;
    }
    
  }

  static void compute_sr_rates(InputSequences &sequences, Partitions &partitions, Tree &tree, std::vector<double> &o_srrates) {
    o_srrates.resize(partitions.size());
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      tree.update_SRcount(partitions[i]);
      partitions[i].normalize_costs(sequences.number() - 1);
      o_srrates[i] = partitions[i].average_sr_rate();
    }
  }

  static void print_stats(const std::string &sequences_file,
                    const std::string &partitions_file,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {
    std::ofstream input_os(output_file.c_str(), std::ofstream::out);
    input_os << "sequences file : " << sequences_file << ::std::endl;
    input_os << "partitions file : " << partitions_file << ::std::endl;
    input_os << "number of tree to generate : " << tree_samples_number << ::std::endl;
    input_os << "number of cpus : " << cpu_number << ::std::endl;
    input_os << "partitions file : " << partitions_file << ::std::endl;

    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences);
    InputPartitions inputPartitions;
    parse_partitions(partitions_file.c_str(), inputPartitions);
    Partitions partitions;
    inputPartitions.generate_partitions(partitions, &sequences);

    input_os << "number of taxa : " << sequences.number() << ::std::endl;
    input_os << "size of sequences : " << sequences.width() << ::std::endl;
    input_os << "number of partitions : " << partitions.size() << ::std::endl;

    Tree tree;
    for (unsigned int i = 0; i < tree_samples_number; ++i) {
      tree.set_random(sequences.number());
      for (unsigned int p = 0; p < partitions.size(); ++p) {
        tree.update_SRcount(partitions[p]);
      }
    }
    for (unsigned int p = 0; p < partitions.size(); ++p) {
      partitions[p].normalize_costs(tree_samples_number * (sequences.number() - 1));
    }
    LoadBalancing lb_naive(partitions, cpu_number);
    LoadBalancing lb_kassian(partitions, cpu_number);
    LoadBalancing lb_weighted(partitions, cpu_number);
    std::cout << "compute naive" << std::endl;
    lb_naive.compute_naive();
    std::cout << "compute kassian" << std::endl;
    lb_kassian.compute_kassian();
    std::cout << "compute weights" << std::endl;
    lb_weighted.compute_kassian_weighted();

//    tree.set_random(sequences.number());
    treatlb(lb_naive, "naive", tree, input_os);    
    treatlb(lb_kassian, "kassian", tree, input_os);    
    treatlb(lb_weighted, "weighted", tree, input_os);    

    
    input_os.close();
  }
  
private:
  static void treatlb(LoadBalancing &lb, const std::string &method, Tree &random_tree, std::ofstream &os) {
    std::cout << "treatlb " << method << std::endl;
    os << "---------------------------------" << std::endl;
    os << method << std::endl;
    os << "Is consistent : " << (lb.is_consistent() ? "yes" : "no") << std::endl;
    os << "Is sites balanced : " << (lb.is_sites_balanced() ? "yes" : "no") << std::endl;
    os << "Is weight balanced : " << (lb.is_weights_balanced() ? "yes" : "no") << std::endl;
    os << "Max partitions difference : " << lb.max_partitions_difference() << std::endl;
    Assignments assignments;
    lb.build_assignments(assignments);
    unsigned int total_sites = 0;
    os << "sites per cpu : ";
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      unsigned int sites = 0;
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        sites += partition.size();
      }
      total_sites += sites;
      os << sites << " ";
    }
    os << "(total : " << total_sites << ")" << std::endl;
    os << "weights per cpu : " ;
    double total_weights = 0;
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      double weights = 0.0;
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        random_tree.update_SRcount(partition);
        partition.normalize_costs(partition.sequences()->number() - 1);
        weights += partition.total_weight();
      }
      total_weights += weights;
      os << weights << " ";
    }
    os << "(total : " << total_weights << ")" << std::endl;
    os << "partitions per cpu : " ; 
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      os << assignments[i].size() << " ";
    }
    os << std::endl;
  }

};

#endif 

