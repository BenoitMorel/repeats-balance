#include "LikelihoodEngine.hpp"  

#include "Tree.hpp"
#include "Partition.hpp"



/* a callback function for performing a full traversal */
static int traverser_full(pll_utree_t * node)
{
  node->data = (void *)~0;
  return 1;
}

LikelihoodEngine::LikelihoodEngine(Tree &tree, Partition &partition):
  tree(tree),
  partition(partition)
{

}
  
void LikelihoodEngine::update_operations()
{
  tree.update_operations(traverser_full);
}
  
void LikelihoodEngine::update_matrices()
{
  partition.update_matrices(tree);
}
  
void LikelihoodEngine::update_partials()
{
  partition.update_partials(tree);
}

double LikelihoodEngine::compute_likelihood()
{
  return partition.compute_likelihood(tree);
}

void LikelihoodEngine::update_sumtable()
{
  partition.update_sumtable(tree);
}

void LikelihoodEngine::compute_derivatives(double *d_f, double *dd_f) 
{	
  partition.compute_derivatives(d_f, dd_f);
}


