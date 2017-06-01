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

bool check_likelihooods(unsigned int states,
    const char *phy_file, 
    const char *part_file,
    std::vector<const char *> newicks,
    std::vector<double> expected_likelihoods,
    bool repeats)
{
  unsigned int attribute = Partition::compute_attribute(repeats, 
		  0, 
		  "avx"); 
  MSA msa(phy_file, states);
  std::vector<Tree *> trees;
  for (unsigned int i = 0; i < newicks.size(); ++i) {
    trees.push_back(new Tree(&msa, newicks[i]));
  }
  LikelihoodEngine engine(trees[0], &msa, part_file, attribute, states, 4, 0); 
  bool result = true;

  for (unsigned int i = 0; i < trees.size(); ++i) {
    engine.set_current_tree(trees[i]);
    engine.update_operations();
    engine.update_matrices();
    engine.update_partials();
    double ll = engine.compute_likelihood();
    if (fabs(ll - expected_likelihoods[i]) > 1.0) {
      std::cerr << "[ERROR] Expected value : " << expected_likelihoods[i] << ", got: " << ll << std::endl;
      result = false;
    }
  }
  for (unsigned int i = 0; i < trees.size(); ++i) {
    delete trees[i];
  }
  return result;
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
  std::vector<const char*> newicks;
  std::vector<double> likelihoods;
  newicks.push_back("data/404_unrooted.newick");
  newicks.push_back("data/404_unrooted9.newick");
  likelihoods.push_back(-379474.0);
  likelihoods.push_back(-379506.0);
  bool ok = true;
  ok &= check_likelihooods(4, "data/404.phy", "data/404.part", newicks, likelihoods, true);
  ok &= check_likelihooods(4, "data/404.phy", "data/404.part", newicks, likelihoods, false);
  return ok;
}

int main(/*int argc, char *params[]*/)
{
  fprintf(stderr, "test_likelihood_404 : %s\n", test_likelihood_404() ? "ok" : "KO");
  fprintf(stderr, "test_likelihood_404_part : %s\n", test_likelihood_404_part() ? "ok" : "KO");
  return 0;

}
