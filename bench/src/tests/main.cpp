#include "../common/repeatsbalance.hpp"
#include <iostream>

double compute_likelihood(unsigned int states,
    const char *newick_file,
    const char *phy_file,
    bool repeats)
{
  unsigned int attribute = Partition::compute_attribute(repeats, 
		  0, 
		  "avx"); 
  Tree tree(newick_file);
  Partition partition(phy_file, tree, attribute, states, 4, 0);
  LikelihoodEngine engine(tree, partition);
  
  engine.update_operations();
  engine.update_matrices();
  engine.update_partials();
  return engine.compute_likelihood();
  
}

bool test_likelihood_404()
{
  double llrepeats = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy", true);
  double lltipinner = compute_likelihood(4, "data/404_unrooted.newick", "data/404.phy", true);
  if (fabs(llrepeats + (double)379474) > 1.0) {
    return false;
  }
  if (fabs(lltipinner + (double)379474) > 1.0) {
    return false;
  }
  return true;
}





int main(/*int argc, char *params[]*/)
{
  fprintf(stderr, "test_likelihood_404 : %s\n", test_likelihood_404() ? "ok" : "KO");
  return 0;

}
