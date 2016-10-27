#ifndef _HELPER_H_
#define _HELPER_H_

#include <string>
#include <fstream>
#include <vector>
#include "repeatsbalance.h"
#include "Partition.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"


struct AssignmentOverview {
  unsigned int cpus;
  std::vector<unsigned int> sites;
  std::vector<unsigned int> partitions;
  std::vector<double> weights;
  unsigned int total_sites;
  double total_weight;
  double max_weight;
  double ratio;
  unsigned int max_partitions;
  unsigned int diff_partitions;

  AssignmentOverview() {
    reset(0);
  }
  void reset(unsigned int cpus) {
    this->cpus = cpus;
    total_sites = max_partitions = diff_partitions = 0;
    total_weight = max_weight = ratio = 0.0;
    sites.resize(cpus);
    std::fill(sites.begin(), sites.end(), 0);
    weights.resize(cpus);
    std::fill(weights.begin(), weights.end(), 0);
    partitions.resize(cpus);
    std::fill(partitions.begin(), partitions.end(), 0);
  }

  friend std::ostream& operator<< (std::ostream &out, const AssignmentOverview &res) {
    out << "Total sites number = " << res.total_sites << std::endl;
    out << "Total weight       = " << res.total_weight << std::endl;
    out << "Max weight         = " << res.max_weight << std::endl;
    out << "Weight ratio       = " << res.ratio << " ((max - average) / max)" << std::endl;
    out << "Max partitions     = " << res.max_partitions << std::endl;
    out << "Sites per cpu      = "; 
    for (unsigned int i = 0; i < res.sites.size(); ++i) { out << res.sites[i] << " ";} out << std::endl;
    out << "Weight per cpu     = "; 
    for (unsigned int i = 0; i < res.weights.size(); ++i) { out << res.weights[i] << " ";} out << std::endl;
    out << "Sites per cpu      = "; 
    for (unsigned int i = 0; i < res.partitions.size(); ++i) { out << res.partitions[i] << " ";} out << std::endl;
    return out;
  }
};

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
    lb_naive.compute_naive();
    lb_kassian.compute_kassian();
    lb_weighted.compute_kassian_weighted();
    
    tree.set_random(sequences.number());

//    tree.set_random(sequences.number());
    AssignmentOverview res_naive, res_kassian, res_weighted;
    treatlb(lb_naive, tree, res_naive);    
    treatlb(lb_kassian, tree, res_kassian);    
    treatlb(lb_weighted, tree, res_weighted);    
    double sites_per_cpu = sequences.width() / cpu_number;
    input_os << std::endl; 
    input_os << "Naive : " << std::endl << res_naive << std::endl;
    input_os << "Kassian : " << std::endl << res_kassian << std::endl;
    input_os << "Kassian weighted : " << std::endl << res_weighted << std::endl;
    

  
    input_os << "Expectations with the sites repeats (if the sites repeats "
                 "computation time is negligible and without the constant "
                 "partition cost)): " << std::endl;
    input_os << "Speedup [no SR, parall] -> // [SR (kassian), parall] : " 
             << sites_per_cpu / res_kassian.max_weight << std::endl;
    input_os << "Speedup [no SR, parall] -> // [SR (kassian weighted), parall] : " 
             << sites_per_cpu / res_weighted.max_weight << std::endl;
    input_os << "Speedup [Kassian, parall] -> [Kassian weighted, parall] : " 
             << res_kassian.max_weight / res_weighted.max_weight << std::endl;
    input_os << "Speedup [no SR, sequential] => [SR (kassian weighted), parall] " 
             << double(res_weighted.total_sites) / res_weighted.max_weight << std::endl;
    input_os.close();
  }
  
private:
  static void treatlb(LoadBalancing &lb, Tree &tree, AssignmentOverview &res) {
    Assignments assignments;
    lb.build_assignments(assignments);
    res.reset(assignments.size());
    res.diff_partitions = lb.max_partitions_difference();
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        tree.update_SRcount(partition);
        partition.normalize_costs(partition.sequences()->number() - 1);
        res.sites[i] += partition.size();
        res.weights[i] += partition.total_weight();
      }
      res.partitions[i] = assign.size();
      res.total_sites += res.sites[i];
      res.total_weight += res.weights[i];
      res.max_weight = std::max(res.weights[i], res.max_weight);
    }
    res.ratio = (res.max_weight - res.total_weight / double(assignments.size())) / res.max_weight;
  }

};

#endif 

