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

    
    std::ofstream input_os((dir + "inputs_infos.txt").c_str(), std::ofstream::out);
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

    LoadBalancing lb_naive(partitions, cpu_number);
    LoadBalancing lb_kassian(partitions, cpu_number);
    lb_naive.compute_naive();
    lb_kassian.compute_kassian();

    tree.set_random(sequences.number());
    treatlb(lb_naive, dir, "naive", tree);    
    treatlb(lb_kassian, dir, "kassian", tree);    

    
    input_os.close();
  }
  
private:
  static void treatlb(LoadBalancing &lb, const std::string &dir, const std::string &method, Tree &random_tree) {
    std::ofstream osdiagnosis((dir + method + "_diagnosis.txt").c_str(), std::ofstream::out);
    std::ofstream ossites((dir + method + "_sitespercpu.txt").c_str(), std::ofstream::out);
    std::ofstream ospartitions((dir + method + "_partitionspercpu.txt").c_str(), std::ofstream::out);
    std::ofstream osaveragePLFC((dir + method + "_averagePLFCpercpu.txt").c_str(), std::ofstream::out);
    osdiagnosis << "Is consistent : " << (lb.is_consistent() ? "yes" : "no") << std::endl;
    osdiagnosis << "Is sites balanced : " << (lb.is_sites_balanced() ? "yes" : "no") << std::endl;
    osdiagnosis << "Max partitions difference : " << lb.max_partitions_difference() << std::endl;
    Assignments assignments;
    lb.build_assignments(assignments);
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      ospartitions << assign.size() << " ";
      unsigned int sites = 0;
      double weights = 0;
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        sites += partition.size();
        random_tree.update_SRcount(partition);
        partition.normalize_costs(partition.sequences()->number());
        weights += partition.total_weight();
      }
      ossites << sites << " ";
      osaveragePLFC << weights << " ";
    }
    ospartitions << std::endl;
    ossites << std::endl;
    osaveragePLFC << std::endl;
  }

};

#endif 

