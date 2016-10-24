#include "LoadBalancing.hpp"
#include "CPU.hpp"
#include "Partition.hpp"
#include <stack>
#include "repeatsbalance.h"
#include <vector>
#include <climits>
#include <cmath>

void LoadBalancing::compute_naive() {
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    _CPUs[i % _CPUs.size()].assign_sites(&_partitions[i], 0, _partitions[i].size());
  }
}

void LoadBalancing::compute_kassian() {
  // Sort the partitions by size in ascending order
  PartitionsPointers sorted_partitions;
  Partition::get_sorted_partitions(_partitions, sorted_partitions);
  // Compute the maximum number of sites per CPU
  unsigned int total_sites_number = 0;
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    total_sites_number += _partitions[i].size();
  }   
  unsigned int max_sites_per_cpu = (total_sites_number - 1) / _CPUs.size() + 1;
  unsigned int r = (max_sites_per_cpu * _CPUs.size()) - total_sites_number;
  unsigned int current_partition_index = 0; // index in sorted_partitons (AND NOT IN _partitions)
  unsigned int full_cpus = 0;
  unsigned int current_cpu = 0;
  // Assign partitions in a cyclic manner to cpus until one is too big 
  for (; current_partition_index < sorted_partitions.size(); ++current_partition_index) {
    const Partition *partition = sorted_partitions[current_partition_index];
    current_cpu = current_partition_index % _CPUs.size();
    if ((partition->size() + _CPUs[current_cpu].weight) > max_sites_per_cpu) {
      // the partition exceeds the current cpu's size, go to the next step of the algo
      break;
    }
    // add the partition !
    _CPUs[current_cpu].assign_sites(partition, 0, partition->size());
    if (_CPUs[current_cpu].weight == max_sites_per_cpu) {
      // one more cpu is exactly full
      if (++full_cpus == (_CPUs.size() - r)) {
        // border case : the remaining cpus should not exceed max_sites_per_cpu - 1
        max_sites_per_cpu--;
      }
      // flag it as full (its harder to rely on max_sites_per_cpu because its value changes)
      _CPUs[current_cpu].tag_full();
    }
  }



  std::stack<CPU *> qlow_;
  std::stack<CPU *> &qlow = qlow_; // hack to assign qhigh to qlow when qlow is empty
  std::stack<CPU *> qhigh;
  for (unsigned int i = 0; i < current_cpu; ++i) {
    if (!_CPUs[current_cpu].is_full()) {
      qhigh.push(&_CPUs[i]);
    }
  }
  for (unsigned int i = current_cpu; i < _CPUs.size(); ++i) {
    if (!_CPUs[current_cpu].is_full()) {
      qlow.push(&_CPUs[i]);
    }
  }
  unsigned int remaining_partition_size = sorted_partitions[current_partition_index]->size();
  while (current_partition_index < _partitions.size() && (qlow.size() || qhigh.size())) {
    const Partition *partition = sorted_partitions[current_partition_index];
    // try to dequeue a process from Qhigh and to fill it 
    if (qhigh.size() && (qhigh.top()->weight + remaining_partition_size >= max_sites_per_cpu)) {
      CPU * cpu = qhigh.top();
      qhigh.pop();
      unsigned int toassign = max_sites_per_cpu - cpu->weight;
      cpu->assign_sites(partition, partition->size() - remaining_partition_size, 
                        toassign);
      remaining_partition_size -= toassign;
      if (++full_cpus == (_CPUs.size() - r)) {
        max_sites_per_cpu--;
      }
    } else if ((qlow.top()->weight + remaining_partition_size >= max_sites_per_cpu)) { // same with qlow
      CPU * cpu = qlow.top();
      qlow.pop();
      unsigned int toassign = max_sites_per_cpu - cpu->weight;
      cpu->assign_sites(partition, partition->size() - remaining_partition_size, 
                        toassign);
      remaining_partition_size -= toassign;
      if (++full_cpus == (_CPUs.size() - r)) {
        max_sites_per_cpu--;
      }
    } else { 
      CPU * cpu = qlow.top();
      qlow.pop();
      cpu->assign_sites(partition, partition->size() - remaining_partition_size,
                        remaining_partition_size);
      remaining_partition_size = 0; 
      qhigh.push(cpu);
    }

    if (!qlow.size()) {
      qlow = qhigh;
    }
    if (!remaining_partition_size) {
      if (++current_partition_index < _partitions.size()) {
        remaining_partition_size = sorted_partitions[current_partition_index]->size();
      }
      
    }
  }
}

void LoadBalancing::compute_kassian_weighted() {

}

bool LoadBalancing::is_consistent() const {
  std::vector<unsigned int> sites_per_partitions(_partitions.size());
  std::fill(sites_per_partitions.begin(), sites_per_partitions.end(), 0);
  for (unsigned int cpu = 0; cpu < _CPUs.size(); ++cpu) {
    for (unsigned int assign = 0; assign < _CPUs[cpu].partitions_number(); ++assign) {
      sites_per_partitions[_CPUs[cpu].partitions[assign]->index()] += _CPUs[cpu].sizes[assign];
    }
  }

  // check that all the correct number of sites is assigned for each partitions
  for (unsigned int partition = 0; partition < _partitions.size(); ++partition) {
    if (sites_per_partitions[partition] != _partitions[partition].size()) {
      SRLOG("Warning : unconsistent load balancing : the number of sites assigned from a partition is invalid");
      return false;
    }
  }

  return true;
}
    
bool LoadBalancing::is_sites_balanced() const {
  unsigned int sites = _CPUs[0].sites_number();
  for (unsigned int i = 1; i < _CPUs.size(); ++i) {
    if (abs(sites - _CPUs[i].sites_number()) > 1) {
      return false;
    }
  } 
  return true;
}


unsigned int LoadBalancing::max_partitions_difference() const {
  unsigned int min = UINT_MAX;
  unsigned int max = 0;
  for (unsigned int i = 0; i < _CPUs.size(); ++i) {
    if (min > _CPUs[i].partitions_number()) {
      min = _CPUs[i].partitions_number();
    }
    if (max < _CPUs[i].partitions_number()) {
      max = _CPUs[i].partitions_number();
    }
  }
  return max - min;
}
