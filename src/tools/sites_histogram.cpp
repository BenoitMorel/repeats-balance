#include "tools.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>
#include <math.h>

typedef std::vector< std::vector<unsigned int> > Histogram;

struct ProcessedPartition {
  MSA *msa;
  Histogram histogram;
  std::vector<char> valid_sites;
};

void compute_histogram(ProcessedPartition &partition)
{
  Histogram &histogram = partition.histogram;
  MSA &msa = *partition.msa;
  histogram.resize(msa.get_length());
  std::fill(histogram.begin(),
      histogram.end(),
      std::vector<unsigned int>(msa.get_states_number() + 1, 0));
  char **seq = msa.get_pll_msa()->sequence;
  for (unsigned int i = 0; i < msa.get_tips_number(); ++i) {
    for (unsigned int j = 0; j < msa.get_length(); ++j) {
      if (seq[i][j] == '-') {
        histogram[j][0]++;
      } else if (seq[i][j] == 'A') {
        histogram[j][1]++;
      } else if (seq[i][j] == 'C') {
        histogram[j][2]++;
      } else if (seq[i][j] == 'G') {
        histogram[j][3]++;
      } else if (seq[i][j] == 'T') {
        histogram[j][4]++;
      }
    }
  }
}

void  display_histogram(Histogram &histogram)
{
  for (unsigned int i = 0; i < histogram.size(); ++i) {
    for (unsigned int j = 0; j < histogram[i].size(); ++j) {
      std::cout << histogram[i][j] << ", ";
    }
    std::cout << std::endl;
  }
}

unsigned int count_unichar(Histogram &histogram, unsigned int lim) {
  unsigned int count = 0;
  for (unsigned int i = 0; i < histogram.size(); ++i) {
    unsigned int count_char = 0;
    for (unsigned int j = 1; j < histogram[i].size(); ++j) {
      count_char += (histogram[i][j] > lim);
    }
    count += (count_char > 1);
  }
  return count;
}

void sites_histogram(int argc, char *params[])
{
  if (argc != 2) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip_file partition_file" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = 4;
  MSA full_msa(phy_filename, states_number);
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(part_filename, initial_partitionning);
  unsigned int partitions_number = initial_partitionning.size();
  std::vector<ProcessedPartition> partitions(partitions_number);
  for (unsigned int i = 0; i < partitions_number; ++i) {
    partitions[i].msa = new MSA(&full_msa, initial_partitionning[i], i);
  }
  
  
  // histograms[i][0] = gaps
  // histograms[i][1] = A
  // histograms[i][2] = C
  // histograms[i][3] = G
  // histograms[i][4] = T 
  for (unsigned int i = 0; i < partitions.size(); ++i) { 
    compute_histogram(partitions[i]);
    display_histogram(partitions[i].histogram);
  }
  
}

