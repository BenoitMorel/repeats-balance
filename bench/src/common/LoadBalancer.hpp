#ifndef LOADBALANCER_HH
#define LOADBALANCER_HH

#include "CoreAssignment.hpp"
#include <stack>
class LoadBalancer {

public:
  void kassian_load_balance(unsigned int cores_number,
      const std::vector<PartitionIntervals> &partitions,
      std::vector<CoreAssignment> &assignments);
  
private:
  void compute_sorted_partitions(const std::vector<PartitionIntervals> &partitions,
    std::vector<const PartitionIntervals *> &sorted_partitions);
  void compute_limit_weight();
  void kassian_first_step(std::vector<CoreAssignment> &assignments);
  void kassian_second_step(std::vector<CoreAssignment> &assignments);
 
  void ksc_initqueues(std::vector<CoreAssignment> &assignments);
  bool can_fill_first_core(std::stack<CoreAssignment*> *stack);
  void fill_first_core(std::stack<CoreAssignment*> *stack);
  void assign_remaining_part();


  unsigned int cores_number;
  std::vector<const PartitionIntervals *> sorted_partitions;
  unsigned int current_part; 
  unsigned int current_core;
  double limit_weight_percore;

  // second step only
  std::stack<CoreAssignment *> *qlow;
  std::stack<CoreAssignment *> *qhigh;
  unsigned curr_part_ass_sites;
};

#endif
