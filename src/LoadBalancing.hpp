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
        CPU &cpu = _CPUs[i % _CPUs.size()];
        cpu.partitions.push_back(&_partitions[i]);
        cpu.offsets.push_back(0);
        cpu.sizes.push_back(_partitions[i].size());
      }
    }

    void compute_kassian() {
       
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
