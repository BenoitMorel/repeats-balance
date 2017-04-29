#include "PartitionIntervals.hpp"
#include <fstream>
#include <string>

void PartitionIntervals::parse(const char *part_file, std::vector<PartitionIntervals> &partitionning)
{
  std::ifstream reader(part_file, std::ifstream::in);
   
  std::string ignore;
  unsigned int start;
  unsigned int end;
  reader >> ignore;
  unsigned int i = 0;
  while(reader) {
    reader >> ignore;
    reader >> ignore;
    reader >> start;
    reader.get();
    reader >> end;
    partitionning.push_back(PartitionIntervals(i++, 1.0));
    partitionning[partitionning.size() - 1].add_interval(start - 1, end - start + 1);
    reader >> ignore;
  }
}

std::ostream& operator<< (std::ostream& os, const PartitionIntervals& partition_intervals)
{
  os << "p" << partition_intervals.get_partition_id() << " ";
  for (unsigned int i = 0; i < partition_intervals.starts.size(); ++i) {
    os << "(" << partition_intervals.starts[i] << ", " << partition_intervals.sizes[i] << ") ";
  }
  return os;
}

void PartitionIntervals::assign_intervals_to(PartitionIntervals &dest,
      unsigned int start, 
      unsigned int size) const
{
  unsigned int remaining = size;
  unsigned int position = 0;
  unsigned int index = 0;
  while (position + sizes[index] <= start) {
    position += sizes[index++];
  }
  unsigned int sites_to_assign = std::min(remaining, (sizes[index] + position) - start);
  dest.add_interval((starts[index] + start) - position, sites_to_assign);
  remaining -= sites_to_assign;
  while (remaining) {
    ++index;
    dest.add_interval(starts[index], std::min(remaining, sizes[index]));
  }
}

