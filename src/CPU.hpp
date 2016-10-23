#ifndef _SB_CPU_H_
#define _SB_CPU_H_

#include <vector> 
#include "Partition.hpp"


struct CPU { 
  std::vector<const Partition *> partitions;
  std::vector<unsigned int> offsets;
  std::vector<unsigned int> sizes;


  void assign_sites(const Partition *partition, 
                      unsigned int offset, 
                      unsigned int size) {
    partitions.push_back(partition);
    offsets.push_back(offset);
    sizes.push_back(size);
  }

  unsigned int partitions_number() const {
    return partitions.size();
  }

  // expansive ! 
  unsigned int sites_number() const {
    int res = 0;
    for (unsigned int i = 0; i < sizes.size(); ++i) {
      res += sizes[i];
    }
    return res;
  }

  // expansive ! 
  unsigned int total_site_costs() const {
    std::cout << "not implemented !" << std::endl;
    return 0;
  }


};

#endif
