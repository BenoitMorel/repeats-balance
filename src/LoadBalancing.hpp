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

    void compute_stupid();

    void compute_kassian();

    void compute_kassian_weighted();

    const std::vector<CPU> &load_balancing() const {
      return _CPUs;
    }

    bool is_consistent() const;

  private:
    const Partitions & _partitions;
    std::vector<CPU> _CPUs;

};

#endif
