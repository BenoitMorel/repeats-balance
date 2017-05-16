#ifndef _LIKELIHOODENGINE_
#define _LIKELIHOODENGINE_

#include "Tree.hpp"
#include "Partition.hpp"
#include "CoreAssignment.hpp"

class LikelihoodEngine {
public:
  
  LikelihoodEngine(const char *newick_file,
    const char *phy_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size);

  LikelihoodEngine(const char *newick_file,
    const char *phy_file,
    const char *part_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size);


  LikelihoodEngine(Tree *tree,
    const std::vector<MSA *> &msas,
    const CoreAssignment &assignment,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size);


  ~LikelihoodEngine();

  void update_operations();

  void update_matrices();

  void update_partials(bool update_repeats = true);

  double compute_likelihood();

  void update_sumtable();

  void compute_derivatives(double *d_f, double *dd_f); 

  Tree &get_tree() {return *tree;}

  std::vector<Partition*> &get_partitions() {return partitions;}

private:
  bool delete_tree;
  Tree *tree;
  std::vector<Partition*> partitions;

};


#endif
