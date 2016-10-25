#include "LoadBalancing.hpp"
#include "Bin.hpp"
#include "Partition.hpp"
#include <stack>
#include "repeatsbalance.h"
#include <vector>
#include <climits>
#include <cmath>

void LoadBalancing::compute_naive() {
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    _bins[i % _bins.size()].assign_sites(&_partitions[i], 0, _partitions[i].size());
  }
}

void LoadBalancing::compute_kassian() {
  // Sort the partitions by size in ascending order
  PartitionsPointers sorted_partitions;
  Partition::get_sorted_partitions(_partitions, sorted_partitions);
  // Compute the maximum number of sites per Bin
  unsigned int total_sites_number = 0;
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    total_sites_number += _partitions[i].size();
  }   
  unsigned int max_sites_per_bin = (total_sites_number - 1) / _bins.size() + 1;
  unsigned int r = (max_sites_per_bin * _bins.size()) - total_sites_number;
  unsigned int current_partition_index = 0; // index in sorted_partitons (AND NOT IN _partitions)
  unsigned int full_bins = 0;
  unsigned int current_bin = 0;
  // Assign partitions in a cyclic manner to bins until one is too big 
  for (; current_partition_index < sorted_partitions.size(); ++current_partition_index) {
    const Partition *partition = sorted_partitions[current_partition_index];
    current_bin = current_partition_index % _bins.size();
    if ((partition->size() + _bins[current_bin].weight) > max_sites_per_bin) {
      // the partition exceeds the current bin's size, go to the next step of the algo
      break;
    }
    // add the partition !
    _bins[current_bin].assign_sites(partition, 0, partition->size());
    if (_bins[current_bin].weight == max_sites_per_bin) {
      // one more bin is exactly full
      if (++full_bins == (_bins.size() - r)) {
        // border case : the remaining bins should not exceed max_sites_per_bin - 1
        max_sites_per_bin--;
      }
      // flag it as full (its harder to rely on max_sites_per_bin because its value changes)
      _bins[current_bin].tag_full();
    }
  }



  std::stack<Bin *> qlow_;
  std::stack<Bin *> &qlow = qlow_; // hack to assign qhigh to qlow when qlow is empty
  std::stack<Bin *> qhigh;
  for (unsigned int i = 0; i < current_bin; ++i) {
    if (!_bins[current_bin].is_full()) {
      qhigh.push(&_bins[i]);
    }
  }
  for (unsigned int i = current_bin; i < _bins.size(); ++i) {
    if (!_bins[current_bin].is_full()) {
      qlow.push(&_bins[i]);
    }
  }
  unsigned int remaining_partition_size = sorted_partitions[current_partition_index]->size();
  while (current_partition_index < _partitions.size() && (qlow.size() || qhigh.size())) {
    const Partition *partition = sorted_partitions[current_partition_index];
    // try to dequeue a process from Qhigh and to fill it 
    if (qhigh.size() && (qhigh.top()->weight + remaining_partition_size >= max_sites_per_bin)) {
      Bin * bin = qhigh.top();
      qhigh.pop();
      unsigned int toassign = max_sites_per_bin - bin->weight;
      bin->assign_sites(partition, partition->size() - remaining_partition_size, 
                        toassign);
      remaining_partition_size -= toassign;
      if (++full_bins == (_bins.size() - r)) {
        max_sites_per_bin--;
      }
    } else if ((qlow.top()->weight + remaining_partition_size >= max_sites_per_bin)) { // same with qlow
      Bin * bin = qlow.top();
      qlow.pop();
      unsigned int toassign = max_sites_per_bin - bin->weight;
      bin->assign_sites(partition, partition->size() - remaining_partition_size, 
                        toassign);
      remaining_partition_size -= toassign;
      if (++full_bins == (_bins.size() - r)) {
        max_sites_per_bin--;
      }
    } else { 
      Bin * bin = qlow.top();
      qlow.pop();
      bin->assign_sites(partition, partition->size() - remaining_partition_size,
                        remaining_partition_size);
      remaining_partition_size = 0; 
      qhigh.push(bin);
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
  const double tolerance = 0.1;
  
  // Sort the partitions by size in ascending order
  PartitionsPointers sorted_partitions;
  Partition::get_sorted_partitions_weighted(_partitions, sorted_partitions);
  // Compute the maximum number of weights per Bin
  double total_weight = 0;
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    total_weight += _partitions[i].total_weight();
  }   
  double max_weight_per_bin = total_weight / _bins.size();
  unsigned int current_partition_index = 0; // index in sorted_partitons (AND NOT IN _partitions)
  unsigned int current_bin = 0;
  // Assign partitions in a cyclic manner to bins until one is too big 
  for (; current_partition_index < sorted_partitions.size(); ++current_partition_index) {
    const Partition *partition = sorted_partitions[current_partition_index];
    current_bin = current_partition_index % _bins.size();
    if ((partition->total_weight() + _bins[current_bin].weight) > max_weight_per_bin) {
      // the partition exceeds the current bin's size, go to the next step of the algo
      break;
    }
    // add the partition !
    _bins[current_bin].assign_sites(partition, 0, partition->size(), partition->total_weight());
    if (_bins[current_bin].weight > max_weight_per_bin + tolerance) {
      _bins[current_bin].tag_full();
    }
  }



  std::stack<Bin *> qlow_;
  std::stack<Bin *> &qlow = qlow_; // hack to assign qhigh to qlow when qlow is empty
  std::stack<Bin *> qhigh;
  for (unsigned int i = 0; i < current_bin; ++i) {
    if (!_bins[current_bin].is_full()) {
      qhigh.push(&_bins[i]);
    }
  }
  for (unsigned int i = current_bin; i < _bins.size(); ++i) {
    if (!_bins[current_bin].is_full()) {
      qlow.push(&_bins[i]);
    }
  }
  double assigned_partition_weight = 0;
  unsigned int assigned_partition_sites = 0;
  while (current_partition_index < _partitions.size() && (qlow.size() || qhigh.size())) {
    const Partition *partition = sorted_partitions[current_partition_index];
    // try to dequeue a process from Qhigh and to fill it 
    if (qhigh.size() && (qhigh.top()->weight + partition->total_weight() - assigned_partition_weight >= max_weight_per_bin)) {
      Bin * bin = qhigh.top();
      qhigh.pop();
      double weight_to_assign = max_weight_per_bin - bin->weight;
      unsigned int sites_toassign = partition->find_upper_bound(weight_to_assign + assigned_partition_weight) - assigned_partition_sites; 
      weight_to_assign = partition->get_accumulated_weight(assigned_partition_sites + sites_toassign) - assigned_partition_weight;
      assigned_partition_weight += weight_to_assign;
      bin->assign_sites(partition, assigned_partition_sites, sites_toassign, weight_to_assign);
      assigned_partition_sites += sites_toassign;
    } else if ((qlow.top()->weight + partition->total_weight() - assigned_partition_weight >= max_weight_per_bin)) { // same with qlow
      Bin * bin = qlow.top();
      qlow.pop();
      double weight_to_assign = max_weight_per_bin - bin->weight;
      unsigned int sites_toassign = partition->find_upper_bound(weight_to_assign + assigned_partition_weight) - assigned_partition_sites; 
      weight_to_assign = partition->get_accumulated_weight(assigned_partition_sites + sites_toassign) - assigned_partition_weight;
      assigned_partition_weight += weight_to_assign;
      bin->assign_sites(partition, assigned_partition_sites, sites_toassign, weight_to_assign);
      assigned_partition_sites += sites_toassign;
    } else { 
      Bin * bin = qlow.top();
      qlow.pop();
      double weight_to_assign = partition->total_weight() - assigned_partition_weight;
      bin->assign_sites(partition, assigned_partition_sites, partition->size() - assigned_partition_sites, weight_to_assign);
      assigned_partition_sites = partition->size();
      qhigh.push(bin);
    }

    if (!qlow.size()) {
      qlow = qhigh;
    }
    if (assigned_partition_sites >= partition->size()) {
      if (++current_partition_index < _partitions.size()) {
        assigned_partition_weight = 0.0;
        assigned_partition_sites = 0;
      }
      
    }
  }
}

void LoadBalancing::build_assignments(Assignments &assignments) {
  assignments = Assignments(_bins.size());
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    Assignment &cpu = assignments[i];
    Bin &bin = _bins[i];
    cpu.resize(bin.partitions_number());
    for (unsigned int j = 0; j < bin.partitions_number(); ++j) {
      const Partition *partition = bin.partitions[j]; 
      cpu[j].init(partition->sequences(), 
                  j,
                  partition->start() + bin.offsets[j],
                  bin.sizes[j]);
    }
  }
}

bool LoadBalancing::is_consistent() const {
  std::vector<unsigned int> sites_per_partitions(_partitions.size());
  std::fill(sites_per_partitions.begin(), sites_per_partitions.end(), 0);
  for (unsigned int bin = 0; bin < _bins.size(); ++bin) {
    for (unsigned int assign = 0; assign < _bins[bin].partitions_number(); ++assign) {
      sites_per_partitions[_bins[bin].partitions[assign]->index()] += _bins[bin].sizes[assign];
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
  unsigned int sites = _bins[0].sites_number();
  for (unsigned int i = 1; i < _bins.size(); ++i) {
    if (abs(sites - _bins[i].sites_number()) > 1) {
      return false;
    }
  } 
  return true;
}

bool LoadBalancing::is_weights_balanced() const {
  double weight = _bins[0].weight;
  for (unsigned int i = 1; i < _bins.size(); ++i) {
    if (fabs(weight - _bins[i].weight) / std::max(weight, _bins[i].weight) > 0.01) {
      std::cout << weight << " " << _bins[i].weight << std::endl;
      return false;
    }
  } 
  return true;
}


unsigned int LoadBalancing::max_partitions_difference() const {
  unsigned int min = UINT_MAX;
  unsigned int max = 0;
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    if (min > _bins[i].partitions_number()) {
      min = _bins[i].partitions_number();
    }
    if (max < _bins[i].partitions_number()) {
      max = _bins[i].partitions_number();
    }
  }
  return max - min;
}
