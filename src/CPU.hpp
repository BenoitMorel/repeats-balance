#ifndef _SB_CPU_H_
#define _SB_CPU_H_

#include <iostream> 
#include <vector> 
#include "Partition.hpp"


struct CPU { 
  std::vector<unsigned int> partitions; // indices of the assigned partitions in the partitions array 
  std::vector<unsigned int> offsets; // offsets off the subpartitions from the start of the partitions
  std::vector<unsigned int> sizes; // sizes off the subpartitions assigned to this cpu


  void assign_sites(unsigned int partition, 
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
