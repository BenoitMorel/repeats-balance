#include "../common/repeatsbalance.hpp"
#include <iostream>

double compute_likelihood(unsigned int states,
    const char *newick_file,
    const char *phy_file,
    const char *part_file,
    bool repeats)
{
  unsigned int attribute = Partition::compute_attribute(repeats, 
		  0, 
		  "avx"); 
  LikelihoodEngine engine(newick_file, phy_file, part_file, attribute, states, 4, 0);
  
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  return engine.compute_likelihood();
  
}

bool check_ll(double expected_ll, double ll) 
{
  if (fabs(expected_ll - ll) > 1.0) {
    std::cerr << "[ERROR] Expected value : " << expected_ll << ", got: " << ll << std::endl;
    return false;
  }
  return true;
}

bool test_likelihood_404()
{
  double llrepeats = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy", 0, true);
  double lltipinner = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy", 0, true);
  return check_ll(-379474.0, llrepeats) && check_ll(-379474.0, lltipinner);
}

bool test_likelihood_404_part()
{
  double llrepeats = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy","data/404.part", true);
  double lltipinner = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy", "data/404.part", true);
  return check_ll(-379474.0, llrepeats) && check_ll(-379474.0, lltipinner);
}


//-379506

int main(/*int argc, char *params[]*/)
{
  fprintf(stderr, "test_likelihood_404 : %s\n", test_likelihood_404() ? "ok" : "KO");
  fprintf(stderr, "test_likelihood_404_part : %s\n", test_likelihood_404_part() ? "ok" : "KO");
  return 0;

}
