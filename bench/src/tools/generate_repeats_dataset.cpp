#include "tools.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <string>      /* printf, scanf, puts, NULL */

class GenMSA {
public:
  GenMSA(unsigned int taxas, const std::string states): 
    states(states),
    labels(taxas),
    msa(taxas) {
    for (unsigned int t = 0; t < taxas_number(); ++t) {
      std::stringstream ss;
      ss << "superfrog" << t;
      labels[t] = ss.str();
    }
  }

  void save(const char *filename) {
    std::ofstream writer(filename);
    writer << taxas_number() << " " << sites_number() << std::endl;
    for (unsigned int t = 0; t < taxas_number(); ++t) {
      writer << labels[t] << " " << msa[t] << std::endl;
    }
  }

  unsigned int taxas_number() const {return msa.size();}
  unsigned int sites_number() const {return msa[0].size();}
public:
    unsigned int taxas;
    std::string states;
    std::vector< std::string > labels;
    std::vector< std::string > msa;
};

class GenPart {
public:
  void add_part(unsigned int length) {
    lengths.push_back(length);
  }

  void save(const char *filename) {
    std::ofstream part_file(filename);
    unsigned int curr_offset = 0;
    for (unsigned int i = 0; i < lengths.size(); ++i) {
      part_file << "DNA, Partition" << i << " = " << curr_offset + 1 << "-" << curr_offset + lengths[i] << std::endl;
      curr_offset += lengths[i];
    }
  }
private:
  std::vector<unsigned int> lengths;
};

void add_repeats_partition(unsigned int sites,
  GenMSA &gen_msa,
  GenPart &gen_part) 
{
  // generate a partition with a lot of repeats 
  std::string & states = gen_msa.states;
  for (unsigned s = 0; s < sites; ++s) {
    unsigned int index = s;
    for (unsigned int t = 0; t < gen_msa.taxas_number(); ++t) {
      gen_msa.msa[t].push_back(states[index % states.size()]);
      index /= states.size();
    }
  }
  gen_part.add_part(sites);
}  
  
void add_random_partition(unsigned int sites,
  GenMSA &gen_msa,
  GenPart &gen_part) 
{
  std::string & states = gen_msa.states;
  for (unsigned int t = 0; t < gen_msa.taxas_number(); ++t) {
    std::string &seq = gen_msa.msa[t];
    for (unsigned s = 0; s < sites; ++s) {
      seq.push_back(states[rand() % states.size()]);
    }
  }
  gen_part.add_part(sites);
}

void generate_repeats_dataset(int argc, char *params[])
{
  if (argc != 3) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "output_phy output_part output_tree" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  //const char *tree_filename = params[i++];

  std::string dna_states; 
  dna_states.push_back('A');
  dna_states.push_back('C');
  dna_states.push_back('G');
  dna_states.push_back('T');
  dna_states.push_back('-');
  unsigned int taxas = 50;
  GenMSA gen_msa(taxas, dna_states);
  GenPart gen_part;

  for (unsigned int i = 0; i < 25; ++i) {
    add_random_partition(10000, gen_msa, gen_part); 
  }
  for (unsigned int i = 0; i < 25; ++i) {
    add_repeats_partition(10000, gen_msa, gen_part); 
  }

  gen_msa.save(phy_filename);
  gen_part.save(part_filename);

}



