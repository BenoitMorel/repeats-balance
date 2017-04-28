#ifndef PARTITION_INTERVALS_HH
#define PARTITION_INTERVALS_HH

/*
 *  Describes the segments of sites that belong to
 *  one partition. Indices start from 0
 */
class PartitionIntervals {
public:
  PartitionIntervals();

  void add_interval(unsigned int start, unsigned int size) 
  {
    starts.push_back(start);
    sizes.push_back(size);
  }

  unsigned int get_intervals_number() const {return starts.size();}
  unsigned int get_start(unsigned int interval_index) const {return starts[interval_index];}
  unsigned int get_size(unsigned int interval_index) const {return sizes[interval_index];}

private:
  std::vector<unsigned int> starts;
  std::vector<unsigned int> sizes;
};


#endif
