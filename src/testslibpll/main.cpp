#include <iostream>
#include "../Helper.hpp"
#include "../TexWriter.hpp"
#include "libpll_helper.h"
#include "pll.h"

#ifndef PLL_ATTRIB_SITES_REPEATS // old pll version
  #define PLL_ATTRIB_SITES_REPEATS 0
#endif

void simulate(const std::string &newick,
              const std::string &partitions_dir,
              unsigned int iterations,
              const std::string &output_tex);

void simulate_trees(const std::string &newicks_path,
                    unsigned int newicks_number,
                    const std::string &partition_dir,
                    unsigned int iterations,
                    const std::string &output_tex); 
void bench(const std::string &newicks_path,
                    unsigned int newicks_number,
                    unsigned int attribute,
                    const std::string &partition_dir,
                    unsigned int iterations);

/**
 *  Compute a load balancing from a given sequence, partitions, cpu numbers, and randomly
 *  generated trees
 *  Export it in an output directory (writes several subsequences files for each cpu)
 * */
void export_lbs() {
    Helper::compute_export_lb("../../data/59/59.phy", "../../data/59/59.part", 30, 10, 
                              "../../results/exports/59_10");        
    Helper::compute_export_lb("../../data/128/128.phy", "../../data/128/128.part", 30, 10, 
                              "../../results/exports/128_10");        
    Helper::compute_export_lb("../../data/404/404.phy", "../../data/404/404.part", 30, 10, 
                              "../../results/exports/404_10");        
    Helper::compute_export_lb("../../data/94/94.phy", "../../data/94/94.part", 30, 10, 
                              "../../results/exports/94_10");        
}

/**
 *  Use the files generated by export_lbs to simulate likelihood computations
 *  on a given tree, amd write the time consumed by each simulated cpu in latex file,
 *  with and without repeats, and with both kassian and weighted methods
 * */
void simulate() {
  simulate("../../data/59/unrooted.newick",
                        "../../results/exports/59_10",
                        100, "../../results/libpll_results/59_10cpus.tex"); 
  simulate("../../data/128/unrooted.newick",
                        "../../results/exports/128_10",
                        100, "../../results/libpll_results/128_10cpus.tex"); 
  simulate("../../data/404/seed42_firsttree.newick",
                        "../../results/exports/404_10",
                        100, "../../results/libpll_results/404_10cpus.tex"); 
}

void simulate_trees() {
  simulate_trees("../../data/404/treechar_samples/tree404.", 2,
                        "../../results/exports/404_10",
                        50, "../../results/libpll_results_pertree/404_10cpus.tex"); 
}
void bench(int argc, char *params[]) {
  if (argc != 8) {
    std::cout << "expected syntaxe : ./test [trees prefix] [partitions dir] [iterations] [trees number] [use repeats] [use tip pattern]i [arch]" << std::endl;
    return;
  }
  std::string data = params[1];
  std::string partitions = params[2];
  unsigned int iterations = atoi(params[3]);
  unsigned int trees_number = atoi(params[4]);
  unsigned int attribute = PLL_ATTRIB_ARCH_CPU;
  if (atoi(params[5])) {
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  }
  if (atoi(params[5]) == 2) {
    attribute |= (1 << 10);;
  }
  if (atoi(params[6])) {
    attribute |= PLL_ATTRIB_PATTERN_TIP;
  }
  if (0 == strcmp("sse", params[7])) {
    attribute |= PLL_ATTRIB_ARCH_SSE;
  }
  if (0 == strcmp("avx", params[7])) {
    attribute |= PLL_ATTRIB_ARCH_AVX;
  } 
  if (0 == strcmp("avx2", params[7])) {
    attribute |= PLL_ATTRIB_ARCH_AVX2;
  }
  if (0 == strcmp("avx2", params[7])) {
    attribute |= PLL_ATTRIB_ARCH_AVX2;
  }
  

  bench(data, trees_number, attribute, partitions, iterations); 
}


int main(int argc, char *argv[]) {
  //export_lbs();
  //simulate();
  //simulate_trees();
  //bench();
  bench(argc, argv);

  return 0;
}


double get_max(const std::vector<double> &v) {
  double res = 0;
  for (unsigned int i = 0; i < v.size(); ++i) {
    res = std::max(res, v[i]);
  }
  return res;
}

void bench(const std::string &newicks_path,
                    unsigned int newicks_number,
                    unsigned int attribute,
                    const std::string &partition_dir,
                    unsigned int iterations) 

{
  for (unsigned int i = 0; i < newicks_number; ++i) {
    std::vector<double> times;
    std::ostringstream os; 
    os << newicks_path << i;    
    compute_loadbalancing(os.str().c_str(), (partition_dir + "/sequential").c_str(), attribute, iterations, times); 
    double sum = 0;
    for (unsigned int j = 0; j < times.size(); ++j) {
      sum += times[j];
    }
    std::cout << sum << " ms" << std::endl;
  }
}


void simulate_trees(const std::string &newicks_path,
                    unsigned int newicks_number,
                    const std::string &partition_dir,
                    unsigned int iterations,
                    const std::string &output_tex) {
  std::vector<std::vector<double> > max_times(4);
  double max = 0;
  double maxsr = 0;
  unsigned attr = PLL_ATTRIB_ARCH_CPU;
  unsigned attr_sr = PLL_ATTRIB_SITES_REPEATS;
  for (unsigned int i = 0; i < newicks_number; ++i) {
    std::vector<std::vector<double> > times(4);
    std::ostringstream os; 
    os << newicks_path << i;    
    compute_loadbalancing(os.str().c_str(), (partition_dir + "/kassian").c_str(), attr, iterations, times[0]); 
    compute_loadbalancing(os.str().c_str(), (partition_dir + "/kassian").c_str(), attr_sr, iterations, times[1]); 
    compute_loadbalancing(os.str().c_str(), (partition_dir + "/weighted").c_str(), attr_sr, iterations, times[2]); 
    compute_loadbalancing(os.str().c_str(), (partition_dir + "/sequential").c_str(), attr_sr, iterations, times[3]); 
    max_times[0].push_back(get_max(times[0]));
    max_times[1].push_back(get_max(times[1]));
    max_times[2].push_back(get_max(times[2]));
    max_times[3].push_back(get_max(times[3]) / times[0].size());
    max = std::max(max, max_times[0][i]);
    maxsr= std::max(maxsr, max_times[1][i]);
    maxsr = std::max(maxsr, max_times[2][i]);
    maxsr = std::max(maxsr, max_times[3][i]);
  }
  max = std::max(max, maxsr);
  TexWriter writer(output_tex);
  writer.write_plot<double>(Plot<double>(maxsr, "Worst cpu time per tree with sites repeats")
                            .add_plot(max_times[1], "green", 'x')
                            .add_plot(max_times[2], "blue", 'x')
                            .add_plot(max_times[3], "red", 'x')
                            .xlabel("Tree")
                            .ylabel("Worst cpu time (ms)")
                            .add_legend("green", "Kassian")
                            .add_legend("red", "Lower bound")
      );
  writer.write_plot<double>(Plot<double>(max, "Libpll default implementation against sites repeats implementation")
                            .add_plot(max_times[0], "green", 'x')
                            .add_plot(max_times[2], "blue", 'x')
                            .xlabel("Tree")
                            .ylabel("Worst cpu time (ms)")
                            .add_legend("green", "Kassian, no sites repeats")
                            .add_legend("blue", "Weighted, sites repeats")
      );
  
}

void simulate(const std::string &newick,
              const std::string &partitions_dir,
              unsigned int iterations,
              const std::string &output_tex) {
  TexWriter writer(output_tex);
  std::vector<double> times;
  double max1, max2;
  unsigned attr = PLL_ATTRIB_ARCH_CPU;
  unsigned attr_sr = PLL_ATTRIB_SITES_REPEATS;
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/kassian").c_str(), attr, iterations, times);
  max1 = get_max(times);
  writer.write_plot<double>(Plot<double>(max1,
                     "Elapsed time to compute the likelihood, without site repeats, with Kassian's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)")); 
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/kassian").c_str(), attr_sr, iterations, times);
  max2 = get_max(times);
  writer.write_plot<double>(Plot<double>(max2,
                     "Elapsed time to compute the likelihood, with site repeats, with Kassian's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)") 
                    .add_horizontal_line(max2, "red", true, "worst cpu time with Kassian")); 
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/weighted").c_str(), attr, iterations, times);
  writer.write_plot<double>(Plot<double>(max1,
                     "Elapsed time to compute the likelihood, without site repeats, with Weighted's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)")
                    );
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/weighted").c_str(), attr_sr, iterations, times);
  writer.write_plot<double>(Plot<double>(max2,
                     "Elapsed time to compute the likelihood, with site repeats, with Weighted's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)") 
                    .add_horizontal_line(max2, "red", true, "worst cpu time with Kassian")); 
  times.clear();
}


