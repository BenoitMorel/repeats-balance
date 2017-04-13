#include "Partition.hpp"
#include <search.h>


const unsigned int  Partition::INVALID_ATTRIBUTE = std::numeric_limits<unsigned int>::max();

Partition::Partition(const char *phy_file, 
  Tree &tree,
  unsigned int attribute_flag, 
  unsigned int states_number,
  unsigned int rate_categories_number,
  unsigned int repeats_lookup_size)
{
  is_dna = (states_number == 4);

  unsigned int tips_number = 0;
  pll_msa_t * msa = pll_phylip_parse_msa(phy_file, &tips_number);

  unsigned int * weight = pll_compress_site_patterns(msa->sequence,
      is_dna ? pll_map_nt : pll_map_aa,
      tips_number,
      &(msa->length));

  partition = pll_partition_create(tips_number,
      tips_number - 2,
      states_number,
      (unsigned int)(msa->length),
      1,
      2 * tips_number - 1,
      rate_categories_number,
      tips_number - 2,
      attribute_flag);
  
  pll_set_pattern_weights(partition, weight);
 
  pll_utree_t *pll_utree = tree.get_pll_tree();
  pll_utree_t ** tips_buffer = (pll_utree_t  **)calloc(tips_number,
      sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(pll_utree, tips_buffer);
  hcreate(tips_number);
  unsigned int * tips_clv_indices = (unsigned int *)malloc(tips_number *
      sizeof(unsigned int));
  for (unsigned int i = 0; i < tips_number; ++i)
  {
    tips_clv_indices[i] = tips_buffer[i]->clv_index;
    ENTRY entry;
    entry.key = tips_buffer[i]->label;
    entry.data = (void *)(tips_clv_indices + i);
    hsearch(entry, ENTER);
  }
  for (unsigned int i = 0; i < tips_number; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    query.key = msa->label[i];
    ENTRY * found = NULL;
    found = hsearch(query,ENTER);
    if (!found) {
      fprintf(stderr, "Sequence with header %s does not appear in the tree\n", msa->label[i]);
    }
    unsigned int tip_clv_index = *((unsigned int *)(found->data));
    pll_set_tip_states(partition, tip_clv_index, 
        is_dna ? pll_map_nt : pll_map_aa,
        msa->sequence[i]);
  }

  if (is_repeats_on() && repeats_lookup_size) {
    set_lookup_size(repeats_lookup_size);  
  }
  
  set_model(Model(states_number, rate_categories_number));

  parameters_indices = std::vector<unsigned int>(rate_categories_number, 0);

  sumtable_buffer = (double *)pll_aligned_alloc(
      partition->sites * partition->rate_cats * partition->states_padded *
      sizeof(double), partition->alignment);
  
  free(weight);
  pll_msa_destroy(msa);
}

Partition::~Partition()
{
  pll_partition_destroy(partition);
  pll_aligned_free(sumtable_buffer);
}

void Partition::set_model(const Model &model)
{
  pll_set_frequencies(partition, 0, &model.frequencies[0]);
  pll_set_subst_params(partition, 0, &model.substitution_parameters[0]);
  pll_set_category_rates(partition, &model.rate_categories[0]); 
}

unsigned int Partition::compute_attribute(unsigned int use_repeats, 
                                          unsigned int additional_attr,
                                          const char *arch)
{
  unsigned int attribute = 0;
  if (!strcmp(arch, "cpu")) {
    attribute = PLL_ATTRIB_ARCH_CPU;
  } else if (!strcmp(arch, "sse")) {
    attribute = PLL_ATTRIB_ARCH_SSE;
  } else if (!strcmp(arch, "avx")) {
    attribute = PLL_ATTRIB_ARCH_AVX;
  } else if (!strcmp(arch, "avx2")) {
    attribute = PLL_ATTRIB_ARCH_AVX2;
  } else {
    return INVALID_ATTRIBUTE;
  }
  if (use_repeats == 1) {
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  } else if(use_repeats == 0) {
    attribute |= PLL_ATTRIB_PATTERN_TIP;
  }
  attribute |= additional_attr;
  return attribute;
}

void Partition::set_lookup_size(unsigned int size)
{
  pll_repeats_t *repeats = partition->repeats;
  repeats->lookup_buffer_size = size;
  free(repeats->lookup_buffer);
  repeats->lookup_buffer = (unsigned int *)
     calloc(repeats->lookup_buffer_size, sizeof(unsigned int));
}

void Partition::update_matrices(const Tree &tree)
{
  pll_update_prob_matrices(partition,
    &parameters_indices[0],
    tree.get_matrix_indices(),
    tree.get_branch_lengths(),
    tree.get_matrix_count());
}

void Partition::update_partials(const Tree &tree)
{
  pll_update_partials(partition, tree.get_operations(), tree.get_operations_number()); 
}

double Partition::compute_likelihood(const Tree &tree) 
{
  const pll_utree_t *pll_tree = tree.get_pll_tree();
  return pll_compute_edge_loglikelihood(partition,
      pll_tree->clv_index,
      pll_tree->scaler_index,
      pll_tree->back->clv_index,
      pll_tree->back->scaler_index,
      pll_tree->pmatrix_index,
      &parameters_indices[0],
      NULL);
}

void Partition::update_sumtable(const Tree &tree)
{
  pll_update_sumtable(partition,
                      tree.get_pll_tree()->clv_index,
                      tree.get_pll_tree()->back->clv_index,
                      PLL_SCALE_BUFFER_NONE,
                      PLL_SCALE_BUFFER_NONE,
                      &parameters_indices[0], 
                      sumtable_buffer);
}

void Partition::compute_derivatives(double *d_f, double *dd_f) 
{
  pll_compute_likelihood_derivatives(partition,
                                     PLL_SCALE_BUFFER_NONE,
                                     PLL_SCALE_BUFFER_NONE,
                                     42.0,
                                     &parameters_indices[0],
                                     sumtable_buffer,
                                     d_f, dd_f);

}



