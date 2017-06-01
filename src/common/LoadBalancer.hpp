#ifndef LOADBALANCER_HH
#define LOADBALANCER_HH

#include "MSA.hpp"
#include "CoreAssignment.hpp"
#include <stack>


struct WeightedMSA {
  const MSA *msa;
  double persite_weight;

  WeightedMSA(): msa(0), persite_weight(1.0) {}
  WeightedMSA(const MSA *msa, double persite_weight): msa(msa), persite_weight(persite_weight) {}

  double get_total_weight() const {return persite_weight * msa->get_length();}
};

class LoadBalancer {

public:
  void compute_weighted_msa(const std::vector<MSA *> &input_msas,
      std::vector<WeightedMSA> &weighted_msa,
      unsigned int pll_attribute);
  
  void kassian_load_balance(unsigned int cores_number,
      const std::vector<WeightedMSA> &input_msas,
      std::vector<CoreAssignment> &assignments);


private:
  void compute_sorted_partitions(const std::vector<WeightedMSA> &partitions,
    std::vector<WeightedMSA> &sorted_partitions);
  void compute_limit_weight();
  void kassian_first_step(std::vector<CoreAssignment> &assignments);
  void kassian_second_step(std::vector<CoreAssignment> &assignments);
 
  void ksc_initqueues(std::vector<CoreAssignment> &assignments);
  bool can_fill_first_core(std::stack<CoreAssignment*> *stack);
  void fill_first_core(std::stack<CoreAssignment*> *stack);
  void assign_remaining_part();


  unsigned int cores_number;
  std::vector<WeightedMSA> sorted_partitions;
  unsigned int current_part; 
  unsigned int current_core;
  double limit_weight_percore;

  // second step only
  std::stack<CoreAssignment *> *qlow;
  std::stack<CoreAssignment *> *qhigh;
  unsigned curr_part_ass_sites;
};

#endif
