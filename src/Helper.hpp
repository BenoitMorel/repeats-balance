#ifndef _HELPER_H_
#define _HELPER_H_

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "repeatsbalance.h"
#include "Partition.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"


struct AssignmentOverview {
  unsigned int cpus;
  std::vector<unsigned int> sites;
  std::vector<unsigned int> partitions;
  std::vector<double> weights;
  unsigned int total_sites;
  unsigned int max_sites;
  double total_weight;
  double max_weight;
  double ratio;
  unsigned int max_partitions;
  unsigned int diff_partitions;
  Assignments assignments;

  AssignmentOverview() {
    reset(0);
  }
  void reset(unsigned int cpus) {
    this->cpus = cpus;
    total_sites = max_partitions = diff_partitions = max_sites = 0;
    total_weight = max_weight = ratio = 0.0;
    sites.resize(cpus);
    std::fill(sites.begin(), sites.end(), 0);
    weights.resize(cpus);
    std::fill(weights.begin(), weights.end(), 0);
    partitions.resize(cpus);
    std::fill(partitions.begin(), partitions.end(), 0);
  }

  friend std::ostream& operator<< (std::ostream &out, const AssignmentOverview &res) {
    out << "Total sites number = " << res.total_sites << std::endl;
    out << "Total weight       = " << res.total_weight << std::endl;
    out << "Max weight         = " << res.max_weight << std::endl;
    out << "Weight ratio       = " << res.ratio << " ((max - average) / max)" << std::endl;
    out << "Max partitions     = " << res.max_partitions << std::endl;
    out << "Sites per cpu      = "; 
    for (unsigned int i = 0; i < res.sites.size(); ++i) { out << res.sites[i] << " ";} out << std::endl;
    out << "Weight per cpu     = "; 
    for (unsigned int i = 0; i < res.weights.size(); ++i) { out << res.weights[i] << " ";} out << std::endl;
    out << "Sites per cpu      = "; 
    for (unsigned int i = 0; i < res.partitions.size(); ++i) { out << res.partitions[i] << " ";} out << std::endl;
    return out;
  }
};

class Helper {
public:

  static void compute_sr_rates(InputSequences &sequences, Partitions &partitions, Tree &tree, std::vector<double> &o_srrates) {
    o_srrates.resize(partitions.size());
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      tree.update_SRcount(partitions[i]);
      partitions[i].normalize_costs(sequences.number() - 1);
      o_srrates[i] = partitions[i].average_sr_rate();
    }
  }
  
  static void experiment4(const std::string &sequences_file,
                    const std::string &partitions_file,
                    const std::string &tree_files,
                    unsigned int trees_number,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {
    treat_trees(sequences_file, partitions_file, tree_files, trees_number, tree_samples_number, cpu_number, output_file);
  }
  
  static void experiment5(const std::string &sequences_file,
                    const std::string &partitions_file,
                    unsigned int trees_number,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {
    treat_trees(sequences_file, partitions_file, "", trees_number, tree_samples_number, cpu_number, output_file);
  }


  static void treat_trees(const std::string &sequences_file,
                    const std::string &partitions_file,
                    const std::string &tree_files,
                    unsigned int trees_number,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {
    std::ofstream os(output_file.c_str());
    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences);
    InputPartitions inputPartitions;
    parse_partitions(partitions_file.c_str(), inputPartitions);
    Partitions partitions;
    inputPartitions.generate_partitions(partitions, &sequences);
    Tree tree_sample;
    for (unsigned int i = 0; i < tree_samples_number; ++i) {
      tree_sample.set_random(sequences.number());
      for (unsigned int p = 0; p < partitions.size(); ++p) {
        tree_sample.update_SRcount(partitions[p]);
      }
    }
    for (unsigned int p = 0; p < partitions.size(); ++p) {
      partitions[p].normalize_costs(tree_samples_number * (sequences.number() - 1));
    }
    LoadBalancing lb_k(partitions, cpu_number);
    LoadBalancing lb(partitions, cpu_number);
    lb_k.compute_kassian();
    lb.compute_kassian_weighted();
    Assignments assignments;
    Assignments assignments_k;
    lb.build_assignments(assignments);
    lb_k.build_assignments(assignments_k);
     
    std::vector<double> average(cpu_number, 0.0);
    std::vector<double> max_weights(trees_number, 0.0);
    std::vector<double> max_weights_k(trees_number, 0.0);
    std::vector<unsigned int> max_weights_indices(trees_number, 0);
    std::vector<std::vector<double> > weights(trees_number, average); // weights[tree][cpu]
    std::vector<std::vector<double> > weights_k(trees_number, average); // weights[tree][cpu]
    std::vector<double> lower_bounds(trees_number, 0.0);
    double lower_bound_average = 0;
    double max = 0;
    for (unsigned int t = 0; t < trees_number; ++t) {
      Tree tree;
      tree.reset();
      if (tree_files.size()) {
        std::ostringstream tmp;
        tmp << tree_files << t;
        parse_tree(tmp.str().c_str(), sequences, tree);
      } else {
          tree.set_random(sequences.number());
      }
      for (unsigned int p = 0; p < partitions.size(); ++p) {
        partitions[p].reset_site_costs();
        tree.update_SRcount(partitions[p]);
        partitions[p].normalize_costs(partitions[p].sequences()->number() - 1);
        lower_bounds[t] += partitions[p].total_weight() / cpu_number;
      }
      lower_bound_average += lower_bounds[t] / trees_number;
      for (unsigned int i = 0; i < assignments_k.size(); ++i) {
        Assignment &assign = assignments_k[i];
        for (unsigned int j = 0; j < assign.size(); ++j) {
          Partition &partition = assign[j];
          partition.reset_site_costs();
          tree.update_SRcount(partition);
          partition.normalize_costs(partition.sequences()->number() - 1);
          weights_k[t][i] += partition.total_weight();
        }
        if (max_weights_k[t] < weights_k[t][i]) {
          max_weights_k[t] = weights_k[t][i];
          max = std::max(max, weights_k[t][i]);
        }
      }
      for (unsigned int i = 0; i < assignments.size(); ++i) {
        Assignment &assign = assignments[i];
        for (unsigned int j = 0; j < assign.size(); ++j) {
          Partition &partition = assign[j];
          partition.reset_site_costs();
          tree.update_SRcount(partition);
          partition.normalize_costs(partition.sequences()->number() - 1);
          weights[t][i] += partition.total_weight();
        }
        average[i] += weights[t][i] / trees_number;
        if (max_weights[t] < weights[t][i]) {
          max_weights[t] = weights[t][i];
          max_weights_indices[t] = i;
          max = std::max(max, weights[t][i]);
        }
      }
    }
    std::vector<double> worse_diff_pertree(trees_number);
    double mean = 0;
    double variance = 0;
    for (unsigned int t = 0; t < trees_number; ++t) {
      worse_diff_pertree[t] = (max_weights[t] - lower_bounds[t]) / lower_bound_average;
      mean += worse_diff_pertree[t] / trees_number;
    }
    for (unsigned int t = 0; t < trees_number; ++t) {
      variance += pow(mean - worse_diff_pertree[t] / trees_number, 2);
    }
    variance = sqrt(variance) / trees_number;
    std::cout << "variance : " << variance << std::endl;
    plot<double>(weights[0], max, 0, lower_bounds[0], 0.0, 0.0, max_weights_indices[0],
      "Weights repartitions on the first tree", (partitions.size() < 50), os);
    plot<double>(weights[trees_number - 1], max, 0, lower_bounds[trees_number - 1], 0.0, 0.0, max_weights_indices[trees_number - 1],
      "Weights repartitions on the last", (partitions.size() < 50), os);
    /*iplot<double>(average, max, 0, lower_bound_average, 0.0, -1, 
      "Average weights repartitions on several trees", (partitions.size() < 50), os);*/
    plot3(max_weights_k, max_weights, lower_bounds, max, os);
    
    /*plot<double>(worse_diff_pertree, (max - lower_bound_average) / lower_bound_average, 0, 0, mean, variance, -1,
      "Per tree, normalized difference between the worse weight and the lower bound", false, os);
     */
  }

  static void experiment3(const std::string &sequences_file,
                          const std::string &output_file) {
    std::ofstream os(output_file.c_str());
    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences); 
    std::vector<double> srvector1(sequences.width());
    std::vector<double> srvector2(sequences.width());
    std::vector<double> acc1(sequences.width());
    std::vector<double> acc2(sequences.width());
    std::fill(srvector1.begin(), srvector1.end(), 0.0);
    std::fill(srvector2.begin(), srvector2.end(), 0.0);
    Tree t1;
    Tree t2;
    t1.set_random(sequences.number());
    t2.set_random(sequences.number());
    t1.update_SRcount(sequences, 0, srvector1);
    t2.update_SRcount(sequences, 0, srvector2);
    std::vector<double> diff(sequences.width());
    // normalize
    for (unsigned int i = 0; i < diff.size(); ++i) {
      srvector1[i] = 1.0 - (srvector1[i] / double(sequences.number()));
      srvector2[i] = 1.0 - (srvector2[i] / double(sequences.number()));
    }
    acc1[0] = srvector1[0];
    acc2[0] = srvector2[0];
    diff[0] = fabs(acc1[0] - acc2[0]);
    for (unsigned int i = 1; i < diff.size(); ++i) {
      acc1[i] = acc1[i-1] + srvector1[i];
      acc2[i] = acc2[i-1] + srvector2[i];
      diff[i] = fabs(acc1[i] - acc2[i]); 
    }
    //plot<double>(srvector1, "srvector computed from t1", os);   
    //plot<double>(srvector2, "srvector computed from t2", os);   
    plot<double>(srvector1, "first PLF-C vector", false, os);   
    plot<double>(srvector2, "second PLF-C vector", false, os);   
    plot<double>(acc1, "one accumulated PLF-C vector", false, os);   
    plot<double>(diff, "difference between two accumulated PLF-C vectors", false, os);   

    os.close();
  }

  /**
  * Test kassian and weighted method with a random independant reference tree
  */
  static void experiment1(const std::string &sequences_file,
                    const std::string &partitions_file,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {
    std::cout << "experiment 1 " << std::endl;
    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences);
    InputPartitions inputPartitions;
    parse_partitions(partitions_file.c_str(), inputPartitions);
    Partitions partitions;
    inputPartitions.generate_partitions(partitions, &sequences);
    Tree tree;
    tree.set_random(sequences.number());
    std::ofstream out((output_file + ".txt").c_str(), std::ofstream::out);
    std::ofstream latex_out((output_file + ".tex").c_str(), std::ofstream::out);
    out << "Experiment 1 !" << std::endl;
    out << "sequences file : " << sequences_file << ::std::endl;
    out << "partitions file : " << partitions_file << ::std::endl;
    print_stats(sequences, partitions, tree, tree_samples_number, cpu_number, out, latex_out);
    latex_out.close();
    out.close();

  }

  /**
  * Test kassian and weighted method with a reference tree parsed in a file
  */
  static void experiment2(const std::string &sequences_file,
                    const std::string &partitions_file,
                    const std::string &tree_file,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    const std::string &output_file) {

    InputSequences sequences;
    parse_sequences(sequences_file.c_str(), sequences);
    InputPartitions inputPartitions;
    parse_partitions(partitions_file.c_str(), inputPartitions);
    Partitions partitions;
    inputPartitions.generate_partitions(partitions, &sequences);
    Tree tree;
    parse_tree(tree_file.c_str(), sequences, tree);
    std::ofstream out((output_file + ".txt").c_str(), std::ofstream::out);
    std::ofstream latex_out((output_file + ".tex").c_str(), std::ofstream::out);
    out << "Experiment 2 !" << std::endl;
    out << "sequences file : " << sequences_file << ::std::endl;
    out << "partitions file : " << partitions_file << ::std::endl;
    out << "tree file : " << partitions_file << ::std::endl;
    print_stats(sequences, partitions, tree, tree_samples_number, cpu_number, out, latex_out);
    out.close();

  }

private:
  static void print_stats(const InputSequences &sequences,
                    Partitions &partitions,
                    Tree &reference_tree,
                    unsigned int tree_samples_number,
                    unsigned int cpu_number, 
                    std::ofstream &out,
                    std::ofstream &latex_out) {
    out << "number of tree to generate : " << tree_samples_number << ::std::endl;
    out << "number of cpus : " << cpu_number << ::std::endl;

    out << "number of taxa : " << sequences.number() << ::std::endl;
    out << "size of sequences : " << sequences.width() << ::std::endl;
    out << "number of partitions : " << partitions.size() << ::std::endl;

    // compute the average PLF-C vector
    Tree tree;
    for (unsigned int i = 0; i < tree_samples_number; ++i) {
      tree.set_random(sequences.number());
      for (unsigned int p = 0; p < partitions.size(); ++p) {
        tree.update_SRcount(partitions[p]);
      }
    }
    for (unsigned int p = 0; p < partitions.size(); ++p) {
      partitions[p].normalize_costs(tree_samples_number * (sequences.number() - 1));
    }

    // generate load balancing with the different algorithms
    LoadBalancing lb_naive(partitions, cpu_number);
    LoadBalancing lb_kassian(partitions, cpu_number);
    LoadBalancing lb_weighted(partitions, cpu_number);
    lb_naive.compute_naive();
    lb_kassian.compute_kassian();
    lb_weighted.compute_kassian_weighted();
    
    // get a summary of the results 
    AssignmentOverview res_naive, res_kassian, res_weighted;
    treatlb(lb_naive, reference_tree, res_naive);    
    treatlb(lb_kassian, reference_tree, res_kassian);    
    treatlb(lb_weighted, reference_tree, res_weighted);    

    // and priiiiiint !
    double sites_per_cpu = sequences.width() / cpu_number;
    out << std::endl; 
    out << "Naive : " << std::endl << res_naive << std::endl;
    out << "Kassian : " << std::endl << res_kassian << std::endl;
    out << "Kassian weighted : " << std::endl << res_weighted << std::endl;
    
    out << "Expectations with the sites repeats (if the sites repeats "
                 "computation time is negligible and without the constant "
                 "partition cost)): " << std::endl;
    out << "Speedup [no SR, parall] -> // [SR (kassian), parall] : " 
             << sites_per_cpu / res_kassian.max_weight << std::endl;
    out << "Speedup [no SR, parall] -> // [SR (kassian weighted), parall] : " 
             << sites_per_cpu / res_weighted.max_weight << std::endl;
    out << "Speedup [Kassian, parall] -> [Kassian weighted, parall] : " 
             << res_kassian.max_weight / res_weighted.max_weight << std::endl;
    out << "Speedup [no SR, sequential] => [SR (kassian weighted), parall] " 
             << double(res_weighted.total_sites) / res_weighted.max_weight << std::endl;
    out << std::endl;
    out << std::endl;
    
    double worse_weight = std::max(res_kassian.max_weight, res_weighted.max_weight);
    double lower_bound = res_naive.total_weight / double(cpu_number); 
    plot<unsigned int>(res_kassian.sites, res_kassian.max_sites, 0, 0, 0, 0.0, -1,
      "Sites repartition with Kassian", (partitions.size() < 50), latex_out);
    plot<double>(res_kassian.weights, worse_weight, res_kassian.max_weight, lower_bound, 0, 0.0, -1,
      "PLF-cost repartition with Kassian", (partitions.size() < 50), latex_out);
    plot<unsigned int>(res_weighted.sites, res_weighted.max_sites, 0, 0, 0, 0.0, -1,
      "Sites repartition with Weighted", (partitions.size() < 50), latex_out);
    plot<double>(res_weighted.weights, worse_weight, 0, 0.0,
      res_kassian.max_weight, lower_bound, -1,
      "PLF-cost repartition with Weighted", (partitions.size() < 50), latex_out);
   
  }
  
  static void treatlb(LoadBalancing &lb, Tree &tree, AssignmentOverview &res) {
    Assignments &assignments = res.assignments;
    lb.build_assignments(assignments);
    res.reset(assignments.size());
    res.diff_partitions = lb.max_partitions_difference();
    for (unsigned int i = 0; i < assignments.size(); ++i) {
      Assignment &assign = assignments[i];
      for (unsigned int j = 0; j < assign.size(); ++j) {
        Partition &partition = assign[j];
        tree.update_SRcount(partition);
        partition.normalize_costs(partition.sequences()->number() - 1);
        res.sites[i] += partition.size();
        res.weights[i] += partition.total_weight();
      }
      res.partitions[i] = assign.size();
      res.total_sites += res.sites[i];
      res.total_weight += res.weights[i];
      res.max_weight = std::max(res.weights[i], res.max_weight);
      res.max_sites = std::max(res.sites[i], res.max_sites);
    }
    res.ratio = (res.max_weight - res.total_weight / double(assignments.size())) / res.max_weight;
  }

  
  static void plot_partitions(const Partitions &partitions, const std::vector<Bin> &bins, const std::string &color_prefix, std::ofstream &os) {
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      os << "\\definecolor{" << color_prefix << i << "}{RGB}{";
      os << rand() % 256 << "," << rand() % 256 << "," << rand() % 256 << "}" << std::endl;
    }
    os << "\\begin{tikzpicture} \n"
          "\\begin{axis}[\n"
          "  title={Partitions repartitions with estimated weights},\n"
          "  ybar stacked, ymin=0,\n"  
          "  bar width=3mm,\n"
          "  symbolic x coords={b0";

    for (unsigned int b = 1; b < bins.size(); ++b) {
      os << ",b" << b;
    }
    
    os << "},\n"
          "  xtick=data,\n"
          "  every node near coord/.style={\n"
          "  } ]" << std::endl;

    std::vector<std::vector<double> > partitions_bins(partitions.size());
    std::vector<double> init(bins.size());
    std::fill(init.begin(), init.end(), 0.0);
    std::fill(partitions_bins.begin(), partitions_bins.end(), init);
    for (unsigned int b = 0; b < bins.size(); ++b) {
      int curr = 0;
      for (unsigned int p = 0; p < bins[b].partitions.size(); ++p) {
        unsigned int partition = bins[b].partitions[p]->index();
        partitions_bins[partition][b] = bins[b].weights[curr++];
      }
    }
    for (unsigned int p = 0; p < partitions.size(); ++p) {
      os << "\\addplot [fill=" << color_prefix << p << "] coordinates {";
      for (unsigned int b = 0; b < bins.size(); ++b) {
        os << "({b" << b << "}," << partitions_bins[p][b] << ") ";
      }
      os << "}; " << std::endl;
    }

    os << "\\end{axis}" << std::endl;
    os << "\\end{tikzpicture}" << std::endl;
  }

  template<typename T>
  static void plot(const std::vector<T> &toplot, const std::string &caption, bool histo, std::ofstream &os) {
    T max = 0;
    for (unsigned int i = 0; i < toplot.size(); ++i) {
      max = std::max(max, toplot[i]);
    }
    plot<T>(toplot, max, (T)0, T(0), T(0), 0.0, -1, caption, histo, os);
  }
  
  static void plot3(const std::vector<double> &kassian, 
                   const std::vector<double> &weighted, 
                   const std::vector<double> &lower_band,
                   double ymax,
                   std::ofstream &os) {
    
    os << "\\begin{minipage}{0.49\\textwidth}" << std::endl;
    os << "\\begin{tikzpicture}[scale=0.75]" << std::endl;
      os << "  \\begin{axis}[ymax=" << ymax * 1.1  << ",ymin=0]" << std::endl;
    os << "    \\addplot[color=green, mark=x] coordinates { ";
    for (unsigned int i = 0; i < kassian.size(); ++i) {
      os << "(" << i << "," << kassian[i] << ") ";
    }
    os << "};" << std::endl;
    os << "    \\addplot[color=blue, mark=x] coordinates { ";
    for (unsigned int i = 0; i < weighted.size(); ++i) {
      os << "(" << i << "," << weighted[i] << ") ";
    }
    os << "};" << std::endl;
    os << "    \\addplot[color=red, mark=x] coordinates { ";
    for (unsigned int i = 0; i < lower_band.size(); ++i) {
      os << "(" << i << "," << lower_band[i] << ") ";
    }
    os << "};" << std::endl;

    os << "\\node[draw=black,thick,rounded corners=2pt,above right=2mm] at (0, 0) {%" << std::endl;
    os << "  \\begin{tabular}{@{}r@{ }l@{}}" << std::endl;
    os << "    \\raisebox{2pt}{\\tikz{\\draw[green] (0,0) -- (5mm,0);}}&kassian\\\\" << std::endl;
    os << "    \\raisebox{2pt}{\\tikz{\\draw[blue] (0,0) -- (5mm,0);}}&kassian weighted\\\\" << std::endl;
    os << "    \\raisebox{2pt}{\\tikz{\\draw[red] (0,0) -- (5mm,0);}}&lower band" << std::endl;
    os << "  \\end{tabular}};" << std::endl;


    os << "\\end{axis}" << std::endl;
    os << "\\end{tikzpicture}" << std::endl;
    os << "\\caption*{ Per tree, maximum weights with kassian and our algorithm, and lower band }" << std::endl;
    os << "\\end{minipage}" << std::endl;
  }

  template<typename T>
  static void plot(const std::vector<T> &toplot, 
                   T ymax, T max_kassian, T lower_bound, 
                   T mean,
                   double variance,
                   int index_max,
                   const std::string &caption, bool histo, std::ofstream &os) {
    os << "\\begin{minipage}{0.49\\textwidth}" << std::endl;
    os << "\\begin{tikzpicture}[scale=0.75]" << std::endl;

    if (histo) {
      os << "  \\begin{axis}[ybar interval, ymax=" << T(double(ymax) * 1.1)  << ",ymin=0, minor y tick num = 3]" << std::endl;
      os << "    \\addplot coordinates { ";
    } else {
      os << "  \\begin{axis}[ymax=" << T(double(ymax) * 1.1)  << ",ymin=0]" << std::endl;
      os << "    \\addplot[color=blue, mark=.] coordinates { ";
    }
    for (unsigned int i = 0; i < toplot.size(); ++i) {
      os << "(" << i << "," << toplot[i] << ") ";
    }
    if (histo) {
      // add one value because latex ignores the last one (i dont know why)
      os << "(" << toplot.size() << ", " << ymax / 2 << ") ";
    }
    os << "};" << std::endl;

    if ((double)max_kassian > 0.01) {
      os << "\\draw [red, dashed] ({rel axis cs:0,0}|-{axis cs:0," << max_kassian 
         << "}) -- ({rel axis cs:1,0}|-{axis cs:" << toplot.size() << "," << max_kassian 
         << "}) node [pos=0.5, above] {Max PLF-C with Kassian};" << std::endl;
    }
    if ((double)lower_bound > 0.01) {
      os << "\\draw [red, dashed] ({rel axis cs:0,0}|-{axis cs:0," << lower_bound 
         << "}) -- ({rel axis cs:1,0}|-{axis cs:" << toplot.size() << "," << lower_bound 
         << "}) node [pos=0.5, below] {Lower bound};" << std::endl;
        if (index_max >= 0) {
          os << "\\draw [red, <->] ({axis cs:" << (double)index_max + 0.5
            << "," << lower_bound << "}) -- ({axis cs:" << (double)index_max + 0.5 << 
            "," << toplot[index_max] << "}) node [pos=0.5, left] { " << toplot[index_max] - lower_bound
            << "};" << std::endl;
        }
    }
    std::cout << "mean : " << mean << std::endl;
    if ((double)mean > 0.00001) {
      os << "\\draw [olive, dashed] ({rel axis cs:0,0}|-{axis cs:0," << mean
         << "}) -- ({rel axis cs:1,0}|-{axis cs:" << toplot.size() << "," << mean 
         << "}) node [pos=0.5, below] {Mean : "<< mean  <<"};" << std::endl;
    }
    if ((double)variance > 0.0000001) {
      os << "\\node[red, draw, scale=1.2] at ({rel axis cs:0.5,0.75}) {Variance: :"
         << variance<< "};" << std::endl;
    }
    os << "  \\end{axis}" << std::endl;
    os << "\\end{tikzpicture}" << std::endl;
    os << "\\caption*{" << caption << "}" << std::endl;
    os << "\\end{minipage}" << std::endl;
  }


};

#endif 

