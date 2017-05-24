#include "PartitionIntervals.hpp"
#include <fstream>
#include <string>

void PartitionIntervals::parse(const char *part_file, std::vector<PartitionIntervals> &partitionning)
{
  std::ifstream reader(part_file, std::ifstream::in);
  if (!reader.is_open()) {
    std::cerr << "[Error] Cannot open " << part_file << std::endl;
    return;
  }
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
    partitionning.push_back(PartitionIntervals(i++));
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

