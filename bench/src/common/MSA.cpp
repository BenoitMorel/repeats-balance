#include "MSA.hpp"

MSA::MSA(const char *phy_filename):
  msa(0),
  weights(0),
  msa_idx(0)
{
  unsigned int tips_number = 0;
  msa = pll_phylip_parse_msa(phy_filename, &tips_number);
}


MSA::MSA(const MSA *original_msa, const PartitionIntervals &intervals, unsigned int states_number, unsigned int msa_index) :
  msa(0),
  weights(0),
  msa_idx(msa_index)
{
  if (original_msa->is_compressed()) {
    std::cerr << "[Warning] MSA::MSA original_msa is already compressed, but weights will be ignored" << std::endl;
  }
  const pll_msa_t *full_msa = original_msa->get_pll_msa(); 
  msa->count = full_msa->count;
  msa->length = intervals.get_total_intervals_size();
  msa->sequence = (char **)calloc(msa->count, sizeof(char *));
  msa->label = (char **)calloc(msa->count, sizeof(char *));
  for (int seqidx = 0; seqidx < msa->count; ++seqidx) {
    // copy label
    if (full_msa->label) {
      unsigned int label_length = strlen(full_msa->label[seqidx]) + 1;
      msa->label[seqidx] = (char *)malloc(label_length * sizeof(char));
      memcpy(msa->label[seqidx], full_msa->label[seqidx], label_length); 
    }
    // create subsequence
    if (!full_msa->sequence[seqidx]) {
      std::cerr << "[Error] Missing sequence " << seqidx << " in create_sub_msa " << std::endl;
      return ;
    }

    unsigned int msa_sequence_offset = 0;
    msa->sequence[seqidx] = (char *)malloc(msa->length * sizeof(char));
    for (unsigned int intidx = 0; intidx < intervals.get_intervals_number(); ++intidx) {
      memcpy(msa->sequence[seqidx] + msa_sequence_offset, 
          full_msa->sequence[seqidx] + intervals.get_start(intidx),
          intervals.get_size(intidx) * sizeof(char));
      msa_sequence_offset += intervals.get_size(intidx);
    }
  }
  weights = pll_compress_site_patterns(msa->sequence,
    (states_number == 4) ? pll_map_nt : pll_map_aa,
    msa->count,
    &(msa->length));
}

MSA::~MSA() 
{
  pll_msa_destroy(msa);
  free(weights);
}


