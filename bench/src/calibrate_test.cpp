
#include "common.h"
#include <iostream>
#include "Calibrator.h"
#include <time.h>

void bench(unsigned int sites,
    unsigned int left_sites,
    unsigned int right_sites,
    unsigned int iterations,
    unsigned int lookup_size,
    unsigned int attribute)
{
  std::vector<char> states;
  pll_operation_t operation;
  states.push_back('A');
  states.push_back('C');
  states.push_back('G');
  states.push_back('T');
  pll_partition_t *partition = CA::Calibrator::create_partition(sites,states,left_sites,right_sites, 
      lookup_size, attribute, operation);
  //if (partition->repeats)
  //  CA::Calibrator::print_stats(partition, &operation);
  std::cout << CA::Calibrator::bench_partials(partition, &operation, iterations) << "ms, ";
  std::cout << CA::Calibrator::bench_likelihood(partition, &operation, iterations) << "ms, ";
  std::cout << CA::Calibrator::bench_sumtable(partition, &operation, iterations) << "ms";
  std::cout << std::endl;
  pll_partition_destroy(partition);
}



void calibrate_test(int argc, char *params[])
{
  srand(42);
  unsigned int repeats =  PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_SITES_REPEATS;
  unsigned int noop =  PLL_ATTRIB_ARCH_AVX;
  unsigned int tippattern= PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_PATTERN_TIP;

  unsigned int sites = 10000;
  unsigned int iterations = 5000;
  std::vector<unsigned int> left_sites;
  std::vector<unsigned int> lookup_sizes;
  left_sites.push_back(8);
  left_sites.push_back(sites/10);
  left_sites.push_back(sites/2);
  lookup_sizes.push_back(200);
  lookup_sizes.push_back(2000000);
  lookup_sizes.push_back(20000000);
  std::cout << "Reference tippattern run: " << std::endl;
  bench(sites, sites, 50, iterations, 0, tippattern);
  std::cout << std::endl;
  for (unsigned int ls = 0; ls < lookup_sizes.size(); ++ls) {
    std::cout << "Lookup size : " << lookup_sizes[ls] << std::endl;
    for (unsigned int i = 0; i < left_sites.size(); i++) {
      std::cout << left_sites[i] << " "; //std::endl;
      bench(sites, left_sites[i], (sites * 2) / 3, iterations, lookup_sizes[ls], repeats);
    }
    std::cout << std::endl;
  }
}

