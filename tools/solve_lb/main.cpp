
#include <vector>

#include <iostream>
using namespace std;


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


void evaluate_costs(unsigned int partitions, 
    unsigned int cores,
    const std::vector< std::vector<unsigned int> > &M,
    const std::vector<double> &lambda,
    double alpha,
    std::vector<double> &costs)
{
  costs = std::vector<double> (cores, 0);
  for (unsigned int p = 0; p < partitions; ++p) {
    for (unsigned int c = 0; c < cores; ++c) {
      costs[c] += (M[c][p] != 0) * alpha + M[c][p] * lambda[p];
    }
  }
  std::cout << "new costs : ";
  for (unsigned int c = 0; c < cores; ++c) {
    std::cout << costs[c] << ", ";
  }
  std::cout << std::endl;

}

unsigned int get_min(const std::vector<double> &costs, std::vector< std::vector<unsigned int> > &M, unsigned int p) {
  unsigned int min = 0;
  for (unsigned int c = 0; c < costs.size(); ++c) {
    if (M[c][p] && costs[min] > costs[c]) {
      min = c;
    }
  }
  return min;
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

int main()
{
  const unsigned int partitions = 11;
  const unsigned int cores = 16;

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

  double mycosts[] = {4990, 4733, 4955, 4463, 4220, 3802, 4286, 4581, 4401, 4831, 4245, 4687, 4227, 4821, 5429, 4685};
  std::vector<double> costs;
  costs.assign(mycosts, mycosts + cores);
  std::vector<double> lambdas;
  double alpha;
  compute_lamdas(partitions, cores, M, costs, alpha, lambdas);
  for (unsigned int p = 0; p < partitions; ++p) {
    std::cout << lambdas[p] << ", ";
  }
  std::cout << std::endl;


  evaluate_costs(partitions, cores, M, lambdas, alpha, costs);
  unsigned int worst_core = get_max(costs);
  // find a partition to reduce (or several)
  for (unsigned int p = 0; p < partitions; ++p) {
    if (M[worst_core][p]) {
      // find the min core that shares this partition
      unsigned int best_core = get_min(costs, M, p);
      // exchange
      unsigned int sites = std::min((unsigned int)10, M[worst_core][p]);
      std::cout << "Exchanging " << sites << " sites between core " << worst_core << " and "
        << best_core << " for partition " << p << std::endl;
      M[worst_core][p] -= sites;
      M[best_core][p] += sites;
    }
  }
  evaluate_costs(partitions, cores, M, lambdas, alpha, costs);
  


}
