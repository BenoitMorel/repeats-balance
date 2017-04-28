#ifndef PARTITION_INTERVALS_HH
#define PARTITION_INTERVALS_HH

#include <vector>
#include <iostream>  
/*
 *  Describes the segments of sites that belong to
 *  one partition. Indices start from 0
 */
class PartitionIntervals {
public:
  PartitionIntervals() : total_size(0) {}

  void add_interval(unsigned int start, unsigned int size) 
  {
    starts.push_back(start);
    sizes.push_back(size);
    total_size += size;
  }

  unsigned int get_intervals_number() const {return starts.size();}
  unsigned int get_total_intervals_size() const {return total_size;}
  unsigned int get_start(unsigned int interval_index) const {return starts[interval_index];}
  unsigned int get_size(unsigned int interval_index) const {return sizes[interval_index];}

  static void parse(const char *part_file, std::vector<PartitionIntervals> &partitionning);

  friend std::ostream& operator<< (std::ostream& stream, const PartitionIntervals& partition_intervals);
private:
  std::vector<unsigned int> starts;
  std::vector<unsigned int> sizes;
  unsigned int total_size;
  
};


#endif
