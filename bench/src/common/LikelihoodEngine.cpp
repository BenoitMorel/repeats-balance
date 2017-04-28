#include "LikelihoodEngine.hpp"  

#include "Tree.hpp"
#include "Partition.hpp"



/* a callback function for performing a full traversal */
static int traverser_full(pll_utree_t * node)
{
  node->data = (void *)~0;
  return 1;
}

LikelihoodEngine::LikelihoodEngine(const char *newick_file,
    const char *phy_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size): tree(newick_file)
{
  partitions.push_back(new Partition(phy_file,
        tree,
        attribute_flag,
        states_number,
        rate_categories_number,
        repeats_lookup_size));
}



LikelihoodEngine::LikelihoodEngine(const char *newick_file,
    const char *phy_file,
    const char *part_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size): tree(newick_file)
{
  if (part_file) {
    std::vector<PartitionIntervals> partition_intervals;
    PartitionIntervals::parse(part_file, partition_intervals);
    unsigned int tips_number = 0;
    pll_msa_t * msa = pll_phylip_parse_msa(phy_file, &tips_number);
    if (!msa) {
      std::cout << "Cannot parse msa " << phy_file << std::endl;
    }
    for (unsigned int i = 0; i < partition_intervals.size(); ++i) {
      partitions.push_back(new Partition(msa,
          partition_intervals[i],
          tree,
          attribute_flag,
          states_number,
          rate_categories_number,
          repeats_lookup_size));
    }
    pll_msa_destroy(msa);
  } else {
    partitions.push_back(new Partition(phy_file,
          tree,
          attribute_flag,
          states_number,
          rate_categories_number,
          repeats_lookup_size));
  }
}


 
LikelihoodEngine::~LikelihoodEngine()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    delete partitions[i];
  }
}


void LikelihoodEngine::update_operations()
{
  tree.update_operations(traverser_full);
}
  
void LikelihoodEngine::update_matrices()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_matrices(tree);
  }
}
  
void LikelihoodEngine::update_partials()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_partials(tree);
  }
}

double LikelihoodEngine::compute_likelihood()
{
  double ll = 0;
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    ll += partitions[i]->compute_likelihood(tree);
  }
  return ll;
}

void LikelihoodEngine::update_sumtable()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_sumtable(tree);
  }
}

void LikelihoodEngine::compute_derivatives(double *d_f, double *dd_f) 
{	
  *d_f = *dd_f = 0;
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    double d = 0;
    double dd = 0;
    partitions[i]->compute_derivatives(&d, &dd);
    *d_f += d;
    *dd_f += dd;
  }
}


