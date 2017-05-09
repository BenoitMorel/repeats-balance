#ifndef _MSA_HPP
#define _MSA_HPP

#include "safepll.h"
#include "PartitionIntervals.hpp"


class MSA {

public:

  MSA(const char *phy_filename);

  MSA(const MSA *original_msa, const PartitionIntervals &intervals, unsigned int states_number,  unsigned int msa_index);

  ~MSA();

  void compress();

  bool is_compressed() const {return weights;}

  const pll_msa_t *get_pll_msa() const {return msa;}

  const unsigned int *get_pattern_weights() const {return weights;}

  unsigned int get_msa_index() const {return msa_idx;}
private:
  pll_msa_t *msa;
  unsigned int *weights;
  unsigned int msa_idx;
};


#endif
