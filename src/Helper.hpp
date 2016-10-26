#ifndef _HELPER_H_
#define _HELPER_H_

#include <string>
#include <fstream>
#include <vector>
#include "repeatsbalance.h"
#include "Partition.hpp"
#include "Tree.hpp"

enum Method {
  naive, kassian, siterepeats
};

class Helper {
public:
  static void count_sr(const std::string &sequences_file,
                      const std::string &partitions_file,
                      const std::string &tree_file) {

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

  static void treat(const std::string &sequences_file,
                    const std::string &partitions_file,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_dir) {
    std::string dir = (output_dir[output_dir.size() - 1] == '/') ? output_dir : (output_dir + "/");
    if (system((std::string("mkdir -p ") + dir).c_str())) {
      std::cerr << "Cannot create directory " << dir << std::endl;
      return;
    }

    
    std::ofstream input_os((dir + "results_infos.txt").c_str(), std::ofstream::out);
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
    os << "sites per cpu : ";
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      unsigned int sites = 0;
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        sites += partition.size();
      }
      os << sites << " ";
    }
    os << std::endl;
    os << "weights per cpu : " ;
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      double weights = 0.0;
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        random_tree.update_SRcount(partition);
        partition.normalize_costs(partition.sequences()->number());
        weights += partition.total_weight();
      }
      os << weights << " ";
    }
    os << std::endl;
    os << "partitions per cpu : " ; 
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      os << assignments[i].size() << " ";
    }
    os << std::endl;
  }

};

#endif 

