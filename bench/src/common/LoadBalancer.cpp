#include "LoadBalancer.hpp"
#include <algorithm>


void LoadBalancer::kassian_load_balance(unsigned int _cores_number,
      const std::vector<PartitionIntervals> &partitions,
      std::vector<CoreAssignment> &assignments)
{
  cores_number = _cores_number;
  for (unsigned int i = 0; i < cores_number; ++i) {
    assignments.push_back(CoreAssignment(i));
  }
  current_part = 0; 
  current_core = 0;

  compute_sorted_partitions(partitions, sorted_partitions);
  compute_limit_weight();
  kassian_first_step(assignments);
  kassian_second_step(assignments);

}

bool LoadBalancer::can_fill_first_core(std::stack<CoreAssignment*> *stack)
{
  if (!stack->size()) {
    return false;
  }
  const PartitionIntervals *partition = sorted_partitions[current_part];
  double predicted_weight = stack->top()->get_core_weight() 
    + partition->get_partition_weight()
    - double(curr_part_ass_sites) * partition->get_persite_weight();
  return predicted_weight >= limit_weight_percore;
}

void LoadBalancer::fill_first_core(std::stack<CoreAssignment *> *stack)
{
  CoreAssignment *core = stack->top();
  stack->pop();
  const PartitionIntervals *partition = sorted_partitions[current_part];
  double weight_to_assign = limit_weight_percore - core->get_core_weight();
  unsigned int sites_to_assign = weight_to_assign / 
    partition->get_persite_weight();
  unsigned int available_sites = partition->get_total_intervals_size() - curr_part_ass_sites;
  // in case we exceed the available number of sites
  sites_to_assign = std::min(sites_to_assign, available_sites);
  PartitionIntervals partition_to_assign(partition->get_partition_id(), 
      partition->get_persite_weight());
  partition->assign_intervals_to(partition_to_assign,
      curr_part_ass_sites,
      sites_to_assign);
  core->add_assignment(partition_to_assign);
  curr_part_ass_sites += sites_to_assign;
}
  
void LoadBalancer::assign_remaining_part()
{
  CoreAssignment *core = qlow->top();
  qlow->pop();
  qhigh->push(core);
  const PartitionIntervals *partition = sorted_partitions[current_part];
  PartitionIntervals partition_to_assign(partition->get_partition_id(),
      partition->get_persite_weight());
  unsigned int sites_to_assign = partition->get_total_intervals_size() - curr_part_ass_sites;
  partition->assign_intervals_to(partition_to_assign,
      curr_part_ass_sites,
      sites_to_assign);
  core->add_assignment(partition_to_assign);
  curr_part_ass_sites += sites_to_assign;
}

void LoadBalancer::kassian_second_step(std::vector<CoreAssignment> &assignments)
{
  
  std::stack<CoreAssignment *> qlow_;
  std::stack<CoreAssignment *> qhigh_;
  qlow = &qlow_; 
  qhigh = &qhigh_; 

  ksc_initqueues(assignments);
  
  curr_part_ass_sites = 0;
  while (current_part < sorted_partitions.size() && (qlow->size() || qhigh->size())) {
    if (can_fill_first_core(qhigh)) {
      fill_first_core(qhigh);
    } else if (can_fill_first_core(qlow)) {
      fill_first_core(qlow);
    } else {

      assign_remaining_part();
    }

    if (!qlow->size()) {
      qlow = qhigh;
    }
    if (curr_part_ass_sites >= 
        sorted_partitions[current_part]->get_total_intervals_size()) {
      curr_part_ass_sites = 0;
      current_part++;
    } 
  }
}

bool compare_partitions(const PartitionIntervals* a, 
    const PartitionIntervals* b) { return (*a < *b); }

void LoadBalancer::compute_sorted_partitions(const std::vector<PartitionIntervals> &partitions,
    std::vector<const PartitionIntervals *> &sorted_partitions)
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    sorted_partitions.push_back(&partitions[i]);
  }
  std::sort(sorted_partitions.begin(), sorted_partitions.end(),
     compare_partitions );
}

void LoadBalancer::compute_limit_weight() 
{
  limit_weight_percore = 0.0;
  for (unsigned int i = 0; i < sorted_partitions.size(); ++i) {
    limit_weight_percore += sorted_partitions[i]->get_partition_weight() / double(cores_number);
  }   
  std::cout << "limit weight: " << limit_weight_percore << std::endl;
}

// Assign partitions in a cyclic manner to cores until one is too big 
void LoadBalancer::kassian_first_step(std::vector<CoreAssignment> &assignments)
{
  for (; current_part < sorted_partitions.size(); ++current_part) {
    const PartitionIntervals *partition = sorted_partitions[current_part];
    current_core = current_part % cores_number;
    if ((partition->get_partition_weight() 
          + assignments[current_core].get_core_weight()) 
          > limit_weight_percore) {
      // the partition exceeds the current bin's size, go to the next step of the algo
      break;
    }
    // add the partition !
    assignments[current_core].add_assignment(*partition);
  }
}
  
void LoadBalancer::ksc_initqueues(std::vector<CoreAssignment> &assignments)
{
  for (unsigned int i = 0; i < current_core; ++i) {
    qhigh->push(&assignments[i]);
  }
  for (unsigned int i = current_core; i < assignments.size(); ++i) {
    qlow->push(&assignments[i]);
  }
}

