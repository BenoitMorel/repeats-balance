#ifndef _SB_BIN_H_
#define _SB_BIN_H_

#include "repeatsbalance.h"
#include <iostream> 
#include <vector> 
#include <climits> 

class Partition;

struct Bin { 
  std::vector<const Partition *> partitions; // partitions assigned to this bin 
  std::vector<unsigned int> offsets; // offsets off the subpartitions from the start of the partitions
  std::vector<unsigned int> sizes; // sizes off the subpartitions assigned to this bin
  double weight;

  Bin() : weight(0.0){
  }

  void assign_sites(const Partition * partition, 
                      unsigned int offset, 
                      unsigned int size) {
    if (!size) {
      return;
    }
    partitions.push_back(partition);
    offsets.push_back(offset);
    sizes.push_back(size);
    weight += size;
    //SRLOG("assign " << partition->index() << " offset : " << offset << " size : " << size << " weight is now " << weight);
  }
  
  void assign_sites(const Partition * partition, 
                      unsigned int offset, 
                      unsigned int size, 
                      double weight) {
    if (!size) {
      return;
    }
    partitions.push_back(partition);
    offsets.push_back(offset);
    sizes.push_back(size);
    this->weight += weight;
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

  void tag_full() {
    weight = UINT_MAX;
  }

  bool is_full() const {
    return UINT_MAX == weight;
  }

};

#endif
