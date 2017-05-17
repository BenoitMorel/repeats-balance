#include "LikelihoodEngine.hpp"  

#include "Tree.hpp"
#include "Partition.hpp"



LikelihoodEngine::LikelihoodEngine(const char *newick_file,
    const char *phy_file,
    const char *part_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size): 
  delete_tree(true)
{
  MSA msa(phy_file, states_number);
  tree = new Tree(&msa, newick_file);
  if (part_file) {
    std::vector<PartitionIntervals> partition_intervals;
    PartitionIntervals::parse(part_file, partition_intervals);
    for (unsigned int i = 0; i < partition_intervals.size(); ++i) {
      MSA submsa(&msa, partition_intervals[i], i);
      submsa.compress();
      partitions.push_back(new Partition(&submsa,
          attribute_flag,
          states_number,
          rate_categories_number,
          repeats_lookup_size));
    }
  } else {
    msa.compress();
    partitions.push_back(new Partition(&msa,
          attribute_flag,
          states_number,
          rate_categories_number,
          repeats_lookup_size));
  }
}

LikelihoodEngine::LikelihoodEngine(Tree *tree,
    const MSA *msa,
    const char *part_file,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size): 
  delete_tree(false),
  tree(tree)
{
  std::vector<PartitionIntervals> partition_intervals;
  PartitionIntervals::parse(part_file, partition_intervals);
  for (unsigned int i = 0; i < partition_intervals.size(); ++i) {
    MSA submsa(msa, partition_intervals[i], i);
    submsa.compress();
    partitions.push_back(new Partition(&submsa,
          attribute_flag,
          states_number,
          rate_categories_number,
          repeats_lookup_size));
    }
}

LikelihoodEngine::LikelihoodEngine(Tree *tree,
    const std::vector<MSA *> &msas,
    const CoreAssignment &assignment,
    unsigned int attribute_flag, 
    unsigned int states_number,
    unsigned int rate_categories_number,
    unsigned int repeats_lookup_size): 
  delete_tree(false),
  tree(tree)
{
  const std::vector<PartitionIntervals> &partition_intervals = assignment.get_assignments();
  for (unsigned int i = 0; i < partition_intervals.size(); ++i) {
    const PartitionIntervals &intervals = partition_intervals[i];
    partitions.push_back(new Partition(msas[intervals.get_partition_id()],
        intervals,
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
  if (delete_tree) {
    delete tree;
  }
}


void LikelihoodEngine::update_operations()
{
  tree->update_operations(Tree::traverser_full);
}
  
void LikelihoodEngine::update_matrices()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_matrices(*tree);
  }
}
  
void LikelihoodEngine::update_partials(bool update_repeats)
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_partials(*tree, update_repeats);
  }
}

double LikelihoodEngine::compute_likelihood()
{
  double ll = 0;
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    ll += partitions[i]->compute_likelihood(*tree);
  }
  return ll;
}

void LikelihoodEngine::update_sumtable()
{
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i]->update_sumtable(*tree);
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

void LikelihoodEngine::set_current_tree(Tree *tree)
{
  if (delete_tree) {
    std::cerr << "[Error] LikelihoodEngine::set_current_tree cannot change "
      << "current tree" << std::endl;
    return;
  }
  this->tree = tree;
}

