#ifndef _SB_CPU_H_
#define _SB_CPU_H_

#include "repeatsbalance.h"
#include <iostream> 
#include <vector> 
#include <climits> 

class Partition;

struct CPU { 
  std::vector<const Partition *> partitions; // partitions assigned to this cpu 
  std::vector<unsigned int> offsets; // offsets off the subpartitions from the start of the partitions
  std::vector<unsigned int> sizes; // sizes off the subpartitions assigned to this cpu
  unsigned int weight;

  CPU() : weight(0){
  }

  void assign_sites(const Partition * partition, 
                      unsigned int offset, 
                      unsigned int size) {
    partitions.push_back(partition);
    offsets.push_back(offset);
    sizes.push_back(size);
    weight += size;
    //SRLOG("assign " << partition->index() << " offset : " << offset << " size : " << size << " weight is now " << weight);
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

  void tag_full() {
    weight = UINT_MAX;
  }

  bool is_full() const {
    return UINT_MAX == weight;
  }

};

#endif
