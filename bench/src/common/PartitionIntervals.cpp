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
  while(reader) {
    reader >> ignore;
    reader >> ignore;
    reader >> start;
    reader.get();
    reader >> end;
    partitionning.push_back(PartitionIntervals());
    partitionning[partitionning.size() - 1].add_interval(start - 1, end - start + 1);
    reader >> ignore;
  }
}

std::ostream& operator<< (std::ostream& os, const PartitionIntervals& partition_intervals)
{
  for (unsigned int i = 0; i < partition_intervals.starts.size(); ++i) {
    os << "(" << partition_intervals.starts[i] << ", " << partition_intervals.sizes[i] << ") ";
  }
  return os;
}


