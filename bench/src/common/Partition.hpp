#ifndef PARTITION_
#define PARTITION_

#include <limits>
#include "Model.hpp"
#include "safepll.h"
#include "Tree.hpp"

class Partition {
public:

  Partition(const char *phy_file, 
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size);
  
  Partition(const char *phy_file, 
    Tree &tree,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size);
  
  virtual ~Partition();

  void set_model(const Model &model);

  void update_matrices(const Tree &tree);

  void update_partials(const Tree &tree);

  double compute_likelihood(const Tree &tree);

  void update_sumtable(const Tree &tree);

  void compute_derivatives(double *d_f, double *dd_f); 

  pll_partition_t *get_partition() {return partition;}

public:
  static unsigned int compute_attribute(unsigned int use_repeats, 
                                        unsigned int additional_attr, 
                                        const char *simd_arch);

  static const unsigned int INVALID_ATTRIBUTE;

private:
  void init_partition(pll_msa_t *msa, 
    const std::vector<unsigned int> &tip_indices,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number = 4,
    unsigned int repeats_lookup_size = 0);
  
  bool is_repeats_on() const {return partition->attributes & PLL_ATTRIB_SITES_REPEATS;}
  void set_lookup_size(unsigned int size);

private:
  pll_partition_t *partition;
  std::vector<unsigned int> parameters_indices;
  double * sumtable_buffer;
  bool is_dna;

};



#endif
