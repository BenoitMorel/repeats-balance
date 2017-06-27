#include "tools.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>
#include <math.h>

typedef std::vector< std::vector<unsigned int> > Histogram;
const unsigned int alphabet_size = 16;


struct ProcessedPartition {
  ProcessedPartition(MSA *msa): msa(msa), valid_sites(msa->get_length(), 1) {
    histogram.resize(msa->get_length());
    std::fill(histogram.begin(),
      histogram.end(),
      std::vector<unsigned int>(alphabet_size , 0));
  }

  ~ProcessedPartition() {
    delete msa;
  }


  void compute_histogram()
  {
    char **seq = msa->get_pll_msa()->sequence;
    for (unsigned int i = 0; i < msa->get_tips_number(); ++i) {
      for (unsigned int j = 0; j < sites_number(); ++j) {
       histogram[j][pll_map_nt[(unsigned int)seq[i][j]]]++;
      }
    }
  }

  void  display_histogram()
  {
    for (unsigned int i = 0; i < sites_number(); ++i) {
      for (unsigned int j = 0; j < histogram[i].size(); ++j) {
        std::cout << histogram[i][j] << ", ";
      }
      std::cout << std::endl;
    }
  }

  void invalidate_uninformative_sites(unsigned int lim) 
  {
    for (unsigned int i = 0; i < sites_number(); ++i) {
      unsigned int count_char = 0;
      // skip gaps (last column in histogram)
      for (unsigned int j = 0; j < histogram[i].size() - 1; ++j) {
        count_char += (histogram[i][j] > lim);
      }
      valid_sites[i] &= (count_char > 1);
    }
  }

  unsigned int sites_number() {
    return msa->get_length();
  }

  void write_valid_seq(unsigned int taxon, std::string &buffer) {
    for (unsigned int i = 0; i < sites_number(); ++i) {
      if (valid_sites[i]) {
        buffer.push_back(msa->get_pll_msa()->sequence[taxon][i]);
      }
    }
  }
  
  unsigned int compute_valid_sites_number() {
    unsigned int res = 0;
    for (unsigned int i = 0; i < valid_sites.size(); ++i) {
      res += valid_sites[i];
    }
    return res;
  }

  MSA *msa;
  Histogram histogram;
  std::vector<char> valid_sites;
};



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
  /////////////////////////
  //  PARAMS
  ///////////////////////
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip_file partition_file output limit" 
      << std::endl;
    return ;
  }
  
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  std::string output(params[i++]);
  unsigned int lim = atoi(params[i++]);
  unsigned int states_number = 4;

  ////////////////////////////
  //  GOGO
  ///////////////////////////
  MSA full_msa(phy_filename, states_number);
  std::vector<PartitionIntervals> initial_partitionning;
  PartitionIntervals::parse(part_filename, initial_partitionning);
  unsigned int partitions_number = initial_partitionning.size();
  std::vector<ProcessedPartition *> partitions;
  for (unsigned int i = 0; i < partitions_number; ++i) {
    partitions.push_back(new ProcessedPartition(new MSA(&full_msa, initial_partitionning[i], i)));
  }
  
  for (unsigned int i = 0; i < partitions.size(); ++i) { 
    partitions[i]->compute_histogram();
    partitions[i]->invalidate_uninformative_sites(lim);
  }
 
  std::ofstream os_phy((output + ".phy").c_str());
  std::ofstream os_part((output + ".part").c_str());
  unsigned int valid_sites_number = 0;
  for (unsigned int p = 0; p < partitions.size(); ++p) {
    if (p) {
      os_part << std::endl;
    }
    os_part << "DNA, PARTITION_" << p << " = " << valid_sites_number + 1;
    valid_sites_number += partitions[p]->compute_valid_sites_number();
    os_part << "-" << valid_sites_number;
  }
  
  os_phy << full_msa.get_tips_number() << " " << valid_sites_number << std::endl;
  for (unsigned int t = 0; t < full_msa.get_tips_number(); ++t) {
    std::string buffer;
    for (unsigned int p = 0; p < partitions.size(); ++p) {
      partitions[p]->write_valid_seq(t, buffer);
    }
    os_phy << full_msa.get_pll_msa()->label[t] << " " << buffer << std::endl;
  }
}

