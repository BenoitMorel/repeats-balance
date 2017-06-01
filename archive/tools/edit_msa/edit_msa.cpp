#include <iostream>
#include <fstream>


void generate_partition(const char *input_phy,
    unsigned int partitions_number,
    const char *output_part)
{
  std::ifstream sequence(input_phy);
  unsigned int species;
  unsigned int sites;
  sequence >> species;
  sequence >> sites;
  std::cout << "sites " << sites << std::endl;
  std::ofstream partitions(output_part);
  unsigned int last_offset = 0;
  for (unsigned int i = 0; i < partitions_number; ++i) {
    partitions << "DNA, p" << i << " =";
    partitions << last_offset + 1 << "-";
    last_offset = ((i + 1)  * sites) / (partitions_number);
    partitions << last_offset;
    if (i != partitions_number - 1) {
      partitions << std::endl;
    }
  }
}

void resize_phy(const char *input_phy,
    unsigned int ratio,
    const char *output_phy) {
  std::ifstream sequence(input_phy);
  unsigned int species;
  unsigned int sites;
  sequence >> species;
  sequence >> sites; 
}


int main() {
  generate_partition("../../data/404/404.phy", 30, "/home/morelbt/github/raxml-ng/bin/data/404/404_30.part");
  generate_partition("../../data/404/404.phy", 15, "/home/morelbt/github/raxml-ng/bin/data/404/404_15.part");
  generate_partition("../../data/404/404.phy", 20, "/home/morelbt/github/raxml-ng/bin/data/404/404_20.part");
  generate_partition("../../data/404/404.phy", 40, "/home/morelbt/github/raxml-ng/bin/data/404/404_40part");
  

  return 1;
}
