#include "MSA.hpp"

MSA::MSA(const char *phy_filename, unsigned int states_number):
  msa(0),
  weights(0),
  states_number(states_number),
  msa_idx(0)
{
  pll_phylip_t * fd = pll_phylip_open(phy_filename, pll_map_phylip);
  if (!fd) {
    std::cerr << "[ERROR] MSA::MSA Cannot open " << phy_filename << std::endl;
  }
  msa = pll_phylip_parse_interleaved(fd);
  if (!msa) {
    std::cerr << "[ERROR] MSA::MSA Cannot parse " << phy_filename << std::endl;
  }
  pll_phylip_close(fd);
}


MSA::MSA(const MSA *original_msa, const PartitionIntervals &intervals, unsigned int msa_index) :
  msa(0),
  weights(0),
  states_number(original_msa->get_states_number()),
  msa_idx(msa_index)
{
  if (original_msa->is_compressed()) {
    weights = (unsigned int *)calloc(intervals.get_total_intervals_size(), sizeof(unsigned int));
  }
  msa = (pll_msa_t*)malloc(sizeof(pll_msa_t));
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
    msa->sequence[seqidx] = (char *)calloc((msa->length + 1), sizeof(char));
    for (unsigned int intidx = 0; intidx < intervals.get_intervals_number(); ++intidx) {
      memcpy(msa->sequence[seqidx] + msa_sequence_offset, 
          full_msa->sequence[seqidx] + intervals.get_start(intidx),
          intervals.get_size(intidx) * sizeof(char));
      // also copy rates 
      if (original_msa->is_compressed() && !seqidx) {
        memcpy(weights + msa_sequence_offset,
            original_msa->weights + intervals.get_start(intidx),
            intervals.get_size(intidx) * sizeof(unsigned int));
      }
      msa_sequence_offset += intervals.get_size(intidx);
    }
  }
}

MSA::~MSA() 
{
  pll_msa_destroy(msa);
  free(weights);
}

void MSA::compress()
{
  if (!msa) {
    std::cerr << "[Error] Null MSA in compress" << std::endl;
  }
  if (is_compressed()) {
    std::cerr << "[Warning] MSA::compress msa is already compressed" << std::endl;
    return;
  }
  weights = pll_compress_site_patterns(msa->sequence,
    (states_number == 4) ? pll_map_nt : pll_map_aa,
    msa->count,
    &(msa->length));
}


