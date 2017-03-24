
#include <vector>

#include <iostream>


void compute_lamdas(unsigned int partitions,
    unsigned int cores,
    std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &costs,
    double &alpha,
    std::vector<double> &lambdas)
{
  lambdas = std::vector<double> (partitions, 0.0);  
  std::vector<unsigned int> percore_partitions(cores, 0);
  std::vector<unsigned int> percore_sites(cores, 0);
  std::vector<unsigned int> perpartition_sites(partitions, 0);
  for (unsigned int p = 0; p < partitions; ++p) {
    for (unsigned int c = 0; c < cores; ++c) {
      percore_partitions[c] += (M[c][p] != 0);
      percore_sites[c] += M[c][p];
      perpartition_sites[p] += M[c][p];
    }
  }

  std::cerr << "warning, this is wrong in the general case" << std::endl;
  std::vector<unsigned int> one_partitions_cores;
  std::vector<unsigned int> two_partitions_cores;
  for (unsigned int c = 0; c < cores; ++c) {
    if (percore_partitions[c] == 1) {
      one_partitions_cores.push_back(c);
    } else if (percore_partitions[c] == 2) {
      two_partitions_cores.push_back(c);
    } else {
      std::cerr << "ERROR UNHANDLED CASE" << std::endl;
    }
  }
  double average_cost_1 = 0;
  double average_cost_2 = 0;
  for (unsigned int i = 0; i < one_partitions_cores.size(); ++i) {
    average_cost_1 += double(costs[one_partitions_cores[i]]) / double(one_partitions_cores.size());
  }
  for (unsigned int i = 0; i < two_partitions_cores.size(); ++i) {
    average_cost_2 += double(costs[two_partitions_cores[i]]) / double(two_partitions_cores.size());
  }

  std::vector<double> corrected_costs(costs);
  alpha = average_cost_2 - average_cost_1;
  for (unsigned int i = 0; i < one_partitions_cores.size(); ++i) {
    corrected_costs[one_partitions_cores[i]] -= alpha;
  }
  for (unsigned int i = 0; i < two_partitions_cores.size(); ++i) {
    corrected_costs[two_partitions_cores[i]] -= 2 * alpha;
  }

  for (unsigned int c = 0; c < cores; ++c) {
    for (unsigned int p = 0; p < partitions; ++p) {
      lambdas[p] += corrected_costs[c] * M[c][p] / (perpartition_sites[p] * percore_sites[c]);
    }
  }
  /*
  for (unsigned int p = 0; p < partitions; ++p) {
    lambdas[p] = pow(lambdas[p], 0.5);
  }
  */
}


void evaluate_costs(const std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &lambda,
    double alpha,
    std::vector<double> &costs)
{
  unsigned int cores = M.size();
  unsigned int partitions = M[0].size();
  costs = std::vector<double> (cores, 0);
  for (unsigned int p = 0; p < partitions; ++p) {
    for (unsigned int c = 0; c < cores; ++c) {
      costs[c] += (M[c][p] != 0) * alpha + M[c][p] * lambda[p];
    }
  }
}

double busy_ratio(const std::vector<double> &costs) {
  double max = 0;
  double sum = 0;
  for (auto c: costs) {
    max = std::max(c, max);
    sum += c;
  }
  return sum / (max * costs.size());
}

  
void print_costs(const std::vector<double> &costs) {
  std::cout << "Costs : ";
  for (unsigned int c = 0; c < costs.size(); ++c) {
    std::cout << costs[c] << ", ";
  }
  std::cout << std::endl;
  std::cout << "busy ratio: "  << busy_ratio(costs);
  std::cout << std::endl;
}

unsigned int get_min(const std::vector<double> &costs, std::vector< std::vector<unsigned int> > &M, unsigned int p) {
  unsigned int min = 0;
  bool assigned = false;
  for (unsigned int c = 0; c < costs.size(); ++c) {
    if (M[c][p] && (!assigned || costs[min] > costs[c])) {
      min = c;
      assigned = true;
    }
  }
  return assigned ? min : (unsigned int)-1;
}

unsigned int get_max(const std::vector<double> &costs) {
  unsigned int max = 0;
  for (unsigned int c = 0; c < costs.size(); ++c) {
    if (costs[max] < costs[c]) {
      max = c;
    }
  }
  return max;
}


void print_code_M(const std::vector< std::vector<unsigned int> > &M) 
{
  for (unsigned int i = 0; i < M.size(); ++i) {
    for (unsigned int j = 0; j < M[i].size(); ++j) {
      if (M[i][j]) {
        std::cout << "M[" << i << "][" << j << "] = " <<M[i][j] << ";";
        std::cout << std::endl;
      }
    }
  }
}


void print_M(const std::vector< std::vector<unsigned int> > &M) 
{
  for (unsigned int i = 0; i < M.size(); ++i) {
    for (unsigned int j = 0; j < M[i].size(); ++j) {
      std::cout << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
}


void random_swap(std::vector< std::vector<unsigned int> > &M, bool nomerge = false) {
  unsigned int cores = M.size();
  unsigned int partitions = M[0].size();
  unsigned int core1 = rand() % cores;
  unsigned int core2 = rand() % cores;
  if (core1 == core2) {
    return;
  }
  std::vector<unsigned int> partitions_c1;
  std::vector<unsigned int> partitions_c2;
  for (unsigned int p = 0; p < partitions; ++p) {
    if (M[core1][p]) {
      partitions_c1.push_back(p);
    }
    if (M[core2][p]) {
      partitions_c2.push_back(p);
    }
  }
  unsigned int p1 = partitions_c1[rand() % partitions_c1.size()];
  unsigned int p2 = partitions_c2[rand() % partitions_c2.size()];
  if (p1 == p2) {
    return;
  }
  if (nomerge) {
    if (M[core1][p2] | M[core2][p1])
      return;
  }
  unsigned int sites1 = M[core1][p1];
  unsigned int sites2 = M[core2][p2];
  M[core1][p2] += sites2;
  M[core2][p1] += sites1;
  M[core1][p1] = 0;
  M[core2][p2] = 0;
}


void optimize_sites(std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &lambdas,
    double alpha)
{
  unsigned int cores = M.size();
  unsigned int partitions = M[0].size();
  std::vector<double> costs;
  evaluate_costs(M, lambdas, alpha, costs);
  unsigned int worst_core;
  double worst_cost;
  for (unsigned int i = 0; i < 1000; ++i) {
    worst_core = get_max(costs);
    worst_cost = costs[worst_core];
    // find a partition to reduce (or several)
    for (unsigned int p = 0; p < partitions; ++p) {
      if (M[worst_core][p]) {
        // find the min core that shares this partition
        unsigned int best_core = get_min(costs, M, p);
        if (best_core != (unsigned int)-1) {
          unsigned int sites = std::min((unsigned int)1, M[worst_core][p]);
          M[worst_core][p] -= sites;
          M[best_core][p] += sites;
          evaluate_costs(M, lambdas, alpha, costs);
        }
      }
    }
  } 
}


void optimize_matrix(std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &lambdas,
    double alpha,
    double &worst_cost,
    unsigned int iterations)
{
  for (unsigned int i = 0; i < iterations; ++i) {
    std::cout << "current wost cost : " << worst_cost << std::endl;
    unsigned int cores = M.size();
    unsigned int partitions = M[0].size();
    std::vector< std::vector<unsigned int> > proposal(M);
    std::vector<double> costs;
    double new_cost;
    random_swap(proposal);
    
    
    //this is just for log
    evaluate_costs(proposal, lambdas, alpha, costs);
    new_cost = costs[get_max(costs)];
    //std::cout << "new cost after swap : " << new_cost << std::endl;
    
    optimize_sites(proposal, lambdas, alpha);
    evaluate_costs(proposal, lambdas, alpha, costs);
    new_cost = costs[get_max(costs)];
    //std::cout << "new cost after optimize : " << new_cost << std::endl;
    if (new_cost < worst_cost) {
      //std::cout << "FOUND A NEW BEST MATRIX!" << std::endl;
      //print_M(proposal);
      worst_cost = new_cost;
      M = proposal;
    }
  }
}

void test_new_assignments(std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &lambdas,
    double alpha,
    double &worst_cost,
    unsigned int iterations)
{
  std::vector<double> costs;
  for (unsigned int i = 0; i < iterations; ++i) {
    std::vector< std::vector<unsigned int> > proposal(M);
    for (unsigned int j = 0; j < 1000; ++j) {
      random_swap(proposal, true);
    }
    print_M(proposal);
    evaluate_costs(proposal, lambdas, alpha, costs);
    double new_cost = costs[get_max(costs)];
    optimize_matrix(proposal, lambdas, alpha, new_cost, 100);
    if (new_cost < worst_cost) {
      M = proposal;
      worst_cost = new_cost;
    }
  }
}

int main()
{
  const unsigned int partitions = 11;
  const unsigned int cores = 16;
  bool use_lambdas = true;
  std::vector< std::vector<unsigned int> > M(cores);
  std::vector<unsigned int> temp(partitions, 0);
  std::fill(M.begin(), M.end(), temp);
  M[0][8] = 99;
  M[0][7] = 366;
  M[1][3] = 109;
  M[1][1] = 356;
  M[2][4] = 339;
  M[2][1] = 126;
  M[3][9] = 431;
  M[3][1] = 34;
  M[4][0] = 465;
  M[5][0] = 465;
  M[6][0] = 465;
  M[7][2] = 410;
  M[7][0] = 55;
  M[8][2] = 465;
  M[9][10] = 49;
  M[9][2] = 416;
  M[10][10] = 465;
  M[11][6] = 61;
  M[11][10] = 404;
  M[12][6] = 465;
  M[13][5] = 266;
  M[13][6] = 199;
  M[14][7] = 13;
  M[14][5] = 452;
  M[15][1] = 144;
  M[15][7] = 321;

  print_M(M);

  //double mycosts[] = {4990, 4733, 4955, 4463, 4220, 3802, 4286, 4581, 4401, 4831, 4245, 4687, 4227, 4821, 5429, 4685};
   double mycosts[] = {19031, 18426, 20058, 18552, 17965, 16237, 17908, 19266, 17997, 18784, 17654, 19122, 18103, 19037, 20567, 18936};


  std::vector<double> costs;
  costs.assign(mycosts, mycosts + cores);
  std::vector<double> lambdas;
  double alpha;
  compute_lamdas(partitions, cores, M, costs, alpha, lambdas);

  if (!use_lambdas) {
    double average_lambda = 0;
    for (unsigned int i = 0; i < lambdas.size(); ++i) {
      average_lambda += lambdas[i]/lambdas.size();
    }
    for (unsigned int i = 0; i < lambdas.size(); ++i) {
      lambdas[i] = average_lambda;
    }
  }


  optimize_sites(M, lambdas, alpha); 
  evaluate_costs(M, lambdas, alpha, costs);
  double worst_cost = costs[get_max(costs)]; 
 
  // in practise, doesn't seem to give better results
   test_new_assignments(M, lambdas, alpha, worst_cost, 100);


  std::cout << "lambdas: ";
  for (unsigned int p = 0; p < partitions; ++p) {
    std::cout << lambdas[p] << ", ";
  }
  std::cout << std::endl;
  std::cout << "alpha :"  << alpha << std::endl;
  print_code_M(M);
  evaluate_costs(M, lambdas, alpha, costs);
  std::cout << "WORST COST : " << costs[get_max(costs)] << std::endl;
  print_costs(costs);
}
