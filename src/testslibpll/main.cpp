#include <iostream>
#include "../Helper.hpp"
#include "../TexWriter.hpp"
#include "libpll_helper.h"

void simulate(const std::string &newick,
              const std::string &partitions_dir,
              unsigned int iterations,
              const std::string &output_tex);

void export_lbs() {
    Helper::compute_export_lb("../../data/59/59.phy", "../../data/59/59.part", 30, 10, 
                              "../../results/exports/59_10");        
    Helper::compute_export_lb("../../data/128/128.phy", "../../data/128/128.part", 30, 10, 
                              "../../results/exports/128_10");        
    Helper::compute_export_lb("../../data/404/404.phy", "../../data/404/404.part", 30, 10, 
                              "../../results/exports/404_10");        
}

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


double get_max(const std::vector<double> &v) {
  double res = 0;
  for (unsigned int i = 0; i < v.size(); ++i) {
    res = std::max(res, v[i]);
  }
  return res;
}

void simulate(const std::string &newick,
              const std::string &partitions_dir,
              unsigned int iterations,
              const std::string &output_tex) {
  TexWriter writer(output_tex);
  std::vector<double> times;
  double max1, max2;
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/kassian").c_str(), 0, iterations, times);
  max1 = get_max(times);
  writer.write_plot<double>(Plot<double>(max1,
                     "Elapsed time to compute the likelihood, without site repeats, with Kassian's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)")); 
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/kassian").c_str(), 1, iterations, times);
  max2 = get_max(times);
  writer.write_plot<double>(Plot<double>(max2,
                     "Elapsed time to compute the likelihood, with site repeats, with Kassian's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)") 
                    .add_horizontal_line(max2, "red", true, "worst cpu time with Kassian")); 
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/weighted").c_str(), 0, iterations, times);
  writer.write_plot<double>(Plot<double>(max1,
                     "Elapsed time to compute the likelihood, without site repeats, with Weighted's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)")
                    );
  times.clear();
  
  compute_loadbalancing(newick.c_str(),(partitions_dir + "/weighted").c_str(), 1, iterations, times);
  writer.write_plot<double>(Plot<double>(max2,
                     "Elapsed time to compute the likelihood, with site repeats, with Weighted's load balancing")
                    .add_plot(times).is_histo(true).xlabel("CPU")
                    .ylabel("Elapsed time (ms)") 
                    .add_horizontal_line(max2, "red", true, "worst cpu time with Kassian")); 
  times.clear();
}


int main() {

  //export_lbs();
  simulate();
  return 0;
}
