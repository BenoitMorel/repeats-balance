#ifndef _RB_LOADBALANCING_H_
#define _RB_LOADBALANCING_H_

#include "CPU.hpp"
#include "Partition.hpp"
#include <vector>



class LoadBalancing {
  public:
    LoadBalancing(const Partitions & partitions, unsigned int CPU_number) : 
      _partitions(partitions),
      _CPUs(CPU_number)
    {
      
    }

    void compute_stupid() {
      for (unsigned int i = 0; i < _partitions.size(); ++i) {
        _CPUs[i % _CPUs.size()].assign_sites(&_partitions[i], 0, _partitions[i].size());
      }
    }

    void compute_kassian() {
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
      unsigned int current_partition = 0;
      unsigned int full_cpus = 0;
      std::vector<unsigned int> cpu_weights(_CPUs.size()); // current number of sites
      std::fill(cpu_weights.begin(),cpu_weights.end(), 0);
      // Assign partitions in a cyclic manner to cpus until one is too big 
      for (; current_partition < sorted_partitions.size(); ++current_partition) {
        const Partition *partition = sorted_partitions[current_partition];
        unsigned int current_cpu = current_partition % _CPUs.size();
        if ((partition->size() + cpu_weights[current_cpu]) > max_sites_per_cpu) {
          // the partition exceeds the current cpu's size, go to the next step of the algo
          break;
        }
        // add the partition !
        _CPUs[current_cpu].assign_sites(partition, 0, partition->size());
        cpu_weights[current_cpu] += partition->size();
        if (cpu_weights[current_cpu] == max_sites_per_cpu) {
          // one more cpu is exactly full
          if (++full_cpus == (_CPUs.size() - r)) {
            // border case : the remaining cpus should not exceed max_sites_per_cpu - 1
            max_sites_per_cpu--;
          }
        }

      }
      // break the remaining partitions to fill the remaining cpus
	    // I use the notations off Kassian's pseudo code
      unsigned int low = _CPUs.size() - 1;	
      unsigned int l = _CPUs.size() - 1;	
      unsigned int high = 0;
      unsigned int h = 0;
      unsigned int add = 0;// TODO why 0 ??
      while (current_partition < _partitions.size()) {
        while (cpu_weights[l] >= max_sites_per_cpu) {
          l = --low;
        }
        while (cpu_weights[h] >= max_sites_per_cpu) {
          h = ++high;
        }
        // try to dequeue a process from Qmore and to fill it 
        if (cpu_weights[h] + add >= max_sites_per_cpu) {
          //_CPUs.assign_sites 
        }

      }
    }

    void compute_kassian_weighted() {

    }

    const std::vector<CPU> &load_balancing() const {
      return _CPUs;
    }

  private:
    const Partitions & _partitions;
    std::vector<CPU> _CPUs;

};

#endif
