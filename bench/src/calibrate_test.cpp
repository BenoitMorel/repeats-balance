
#include "common.h"
#include <iostream>
#include "Calibrator.h"
#include <time.h>

void bench(unsigned int sites,
    unsigned int left_sites,
    unsigned int right_sites,
    unsigned int iterations,
    unsigned int attribute)
{
  std::vector<char> states;
  pll_operation_t operation;
  states.push_back('A');
  states.push_back('C');
  states.push_back('G');
  states.push_back('T');
  pll_partition_t *partition = CA::Calibrator::create_partition(sites,states,left_sites,right_sites, 
      attribute, operation);
  //if (partition->repeats)
  //  CA::Calibrator::print_stats(partition, &operation);
  std::cout << CA::Calibrator::bench_partials(partition, &operation, iterations) << "ms, ";
  //std::cout << CA::Calibrator::bench_likelihood(partition, &operation, iterations) << "ms, ";
  //std::cout << CA::Calibrator::bench_sumtable(partition, &operation, iterations) << "ms";
  std::cout << std::endl;
  pll_partition_destroy(partition);
}



void calibrate_test(int argc, char *params[])
{
  srand(42);
  unsigned int repeats =  PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_SITES_REPEATS;
  unsigned int noop =  PLL_ATTRIB_ARCH_AVX;
  unsigned int tippattern= PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_PATTERN_TIP;

  unsigned int sites = 7429;
  unsigned int iterations = 10000;
  std::vector<unsigned int> left_sites;
  left_sites.push_back(8);
  left_sites.push_back(32);
  left_sites.push_back(50);
  left_sites.push_back(sites/80);
  left_sites.push_back(sites/40);
  left_sites.push_back(sites/20);
  left_sites.push_back(sites/10);
  left_sites.push_back(sites/5);
  left_sites.push_back(sites/2);
  left_sites.push_back(sites);
  bench(sites, sites, 50, iterations, tippattern);
  for (unsigned int i = 0; i < left_sites.size(); i++) {
    std::cout << left_sites[i] << " "; //std::endl;
    bench(sites, sites, left_sites[i], iterations, repeats);
  }
}

