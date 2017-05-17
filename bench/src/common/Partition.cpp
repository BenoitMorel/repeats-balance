#include "Partition.hpp"


const unsigned int  Partition::INVALID_ATTRIBUTE = std::numeric_limits<unsigned int>::max();


Partition::Partition(const MSA *compressed_msa, 
  unsigned int attribute_flag, 
  unsigned int states_number,
  unsigned int rate_categories_number,
  unsigned int repeats_lookup_size)
{
  init_partition(compressed_msa->get_pll_msa(), compressed_msa->get_pattern_weights(),  attribute_flag, states_number,
    rate_categories_number, repeats_lookup_size);
}

Partition::Partition(const MSA *compressed_msa, 
  const PartitionIntervals &intervals,
  unsigned int attribute_flag, 
  unsigned int states_number,
  unsigned int rate_categories_number,
  unsigned int repeats_lookup_size)
{
  MSA submsa(compressed_msa, intervals, compressed_msa->get_msa_index());
  init_partition(submsa.get_pll_msa(), submsa.get_pattern_weights(), attribute_flag, states_number,
    rate_categories_number, repeats_lookup_size);
}
  

void Partition::init_partition(const pll_msa_t *compressed_msa, 
  const unsigned int *weights,
  unsigned int attribute_flag, 
  unsigned int states_number,
  unsigned int rate_categories_number,
  unsigned int repeats_lookup_size)
{
  is_dna = (states_number == 4);

  unsigned int tips_number = compressed_msa->count;

  partition = pll_partition_create(tips_number,
      tips_number - 2,
      states_number,
      (unsigned int)(compressed_msa->length),
      1,
      2 * tips_number - 1,
      rate_categories_number,
      tips_number - 2,
      attribute_flag);
  
  pll_set_pattern_weights(partition, weights);
  for (unsigned int i = 0; i < tips_number; ++i) {
    pll_set_tip_states(partition, i, 
      is_dna ? pll_map_nt : pll_map_aa,
      compressed_msa->sequence[i]);
  }

  if (is_repeats_on() && repeats_lookup_size) {
    set_lookup_size(repeats_lookup_size);  
  }
  
  set_model(Model(states_number, rate_categories_number));

  parameters_indices = std::vector<unsigned int>(rate_categories_number, 0);

  sumtable_buffer = (double *)pll_aligned_alloc(
      partition->sites * partition->rate_cats * partition->states_padded *
      sizeof(double), partition->alignment);
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

void Partition::update_partials(const Tree &tree, bool update_repeats)
{
  pll_update_partials_rep(partition, tree.get_operations(), tree.get_operations_number(), update_repeats); 
}

double Partition::compute_likelihood(const Tree &tree) 
{
  const pll_unode_t *pll_root = tree.get_pll_root();
  return pll_compute_edge_loglikelihood(partition,
      pll_root->clv_index,
      pll_root->scaler_index,
      pll_root->back->clv_index,
      pll_root->back->scaler_index,
      pll_root->pmatrix_index,
      &parameters_indices[0],
      NULL);
}

void Partition::update_sumtable(const Tree &tree)
{
  pll_update_sumtable(partition,
                      tree.get_pll_root()->clv_index,
                      tree.get_pll_root()->back->clv_index,
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

void Partition::create_sub_msa(const pll_msa_t *msa, const unsigned int *weights, const PartitionIntervals &intervals, 
      pll_msa_t *&submsa, unsigned int *&subweights) 
{
  submsa = 0;
  subweights = 0;
  // usee malloc because pll_msa_destroy uses free
  if (!msa || !msa->sequence) {
    std::cerr << "[Error] No msa sequence" << std::endl;
    return;
  }
  submsa = (pll_msa_t *)malloc(sizeof(pll_msa_t));  
  submsa->count = msa->count;
  submsa->length = intervals.get_total_intervals_size();
  submsa->sequence = (char **)calloc(submsa->count, sizeof(char *));
  submsa->label = (char **)calloc(submsa->count, sizeof(char *));
  if (weights) {
    subweights = (unsigned int*)calloc(submsa->length, sizeof(unsigned int));
  }
  for (int seqidx = 0; seqidx < submsa->count; ++seqidx) {
    // copy label
    if (msa->label) {
      unsigned int label_length = strlen(msa->label[seqidx]) + 1;
      submsa->label[seqidx] = (char *)malloc(label_length * sizeof(char));
      memcpy(submsa->label[seqidx], msa->label[seqidx], label_length); 
    }
    // create subsequence
    if (!msa->sequence[seqidx]) {
      std::cerr << "[Error] Missing sequence " << seqidx << " in create_sub_msa " << std::endl;
      return ;
    }

    unsigned int submsa_sequence_offset = 0;
    submsa->sequence[seqidx] = (char *)malloc(submsa->length * sizeof(char));
    for (unsigned int intidx = 0; intidx < intervals.get_intervals_number(); ++intidx) {
      memcpy(submsa->sequence[seqidx] + submsa_sequence_offset, 
          msa->sequence[seqidx] + intervals.get_start(intidx),
          intervals.get_size(intidx) * sizeof(char));
      // also copy rates 
      if (!seqidx && weights) {
        memcpy(subweights + submsa_sequence_offset,
            weights + intervals.get_start(intidx),
            intervals.get_size(intidx) * sizeof(unsigned int));
      }
      submsa_sequence_offset += intervals.get_size(intidx);
    }
  }
}


double Partition::get_unique_repeats_pattern_ratio() const
{
  if (!partition->repeats) {
    return 1.0;
  }
  unsigned int total_patterns = 0;
  for (unsigned int node = partition->tips; node < partition->tips + partition->clv_buffers; ++node) {
    unsigned int node_patterns = partition->repeats->pernode_max_id[node];
    node_patterns = node_patterns ? node_patterns : partition->sites;
    total_patterns += node_patterns;
  }
  return double(total_patterns) / double(partition->sites * partition->clv_buffers);
}

