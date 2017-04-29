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
  PartitionIntervals(unsigned int partition_id,
      double persite_weight) :
    total_size(0),
    partition_id(partition_id),
    persite_weight(persite_weight) {}

  void add_interval(unsigned int start, unsigned int size) 
  {
    starts.push_back(start);
    sizes.push_back(size);
    total_size += size;
  }

  unsigned int get_partition_id() const {return partition_id;}
  unsigned int get_intervals_number() const {return starts.size();}
  unsigned int get_total_intervals_size() const {return total_size;}
  unsigned int get_start(unsigned int interval_index) const {return starts[interval_index];}
  unsigned int get_size(unsigned int interval_index) const {return sizes[interval_index];}


  // assign size sites from start from this partition
  // to dest
  // warning: implem could be more efficient
  void assign_intervals_to(PartitionIntervals &dest,
      unsigned int start, 
      unsigned int size) const;


  friend std::ostream& operator<< (std::ostream& stream, const PartitionIntervals& partition_intervals);

  double get_partition_weight() const {return persite_weight * double(total_size);}
  double get_persite_weight() const {return persite_weight;}

  bool operator < (const PartitionIntervals& p) const {
    return  get_partition_weight() < p.get_partition_weight();
  }

  static void parse(const char *part_file, std::vector<PartitionIntervals> &partitionning);

private:
  std::vector<unsigned int> starts;
  std::vector<unsigned int> sizes;
  unsigned int total_size;
  unsigned int partition_id;
  double persite_weight; 
};


#endif
