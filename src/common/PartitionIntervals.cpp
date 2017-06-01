#include "PartitionIntervals.hpp"
#include <fstream>
#include <string>

// merge is the number of partitions that should be merged together (usually 1)
unsigned int get_merge(std::ifstream &reader)
{
  unsigned int merge = 1;
  std::string ignore;
  std::string magic;
  reader >> magic;
  reader.seekg(0, std::ios::beg); 
  if (!magic.compare("div")) {
    reader >> ignore;
    reader >> merge;
  }
  return merge;
}


void read_line(std::ifstream &reader, PartitionIntervals &partition_intervals)
{
  std::string ignore;
  unsigned int start;
  unsigned int end;
  reader >> ignore;
  reader >> ignore;
  reader >> ignore;
  
  do {
    reader >> start;
    if ('-' != reader.peek()) {
      partition_intervals.add_interval(start - 1, 1);
      continue;
    }
    reader.get();
    reader >> end;
    if (reader.peek() == '\\') {
      reader.get();
      unsigned int division = 0;
      reader >> division;
      for (unsigned int i = start; i <= end; i += division) {
        partition_intervals.add_interval(i - 1, 1);
      }
    } else {
      partition_intervals.add_interval(start - 1, end - start + 1);
    }
  } while (reader && reader.get() == ','); 
  // unvalidate reader if we are at the end
  reader.get();
  reader.seekg(-2, std::ios::cur);
}

void PartitionIntervals::parse(const char *part_file, std::vector<PartitionIntervals> &partitionning)
{
  std::ifstream reader(part_file, std::ifstream::in);
  if (!reader.is_open()) {
    std::cerr << "[Error] Cannot open " << part_file << std::endl;
    return;
  }
  
  unsigned int merge = get_merge(reader);
  unsigned int i = 0;
  while(reader) {
    if (i++ % merge == 0) {
      partitionning.push_back(PartitionIntervals(i));
    }
    read_line(reader, partitionning[partitionning.size() - 1]);
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

