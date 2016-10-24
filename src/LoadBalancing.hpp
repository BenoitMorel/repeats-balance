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

    void compute_naive();

    void compute_kassian();

    void compute_kassian_weighted();

    const std::vector<CPU> &load_balancing() const {
      return _CPUs;
    }


    // the follozing mehtods are only used to test and make statistics
    bool is_consistent() const;
    bool is_sites_balanced() const;
    unsigned int max_partitions_difference() const; // different between the max and the min number of partitions in the cpus
  private:
    const Partitions & _partitions;
    std::vector<CPU> _CPUs;

};

#endif
