#include "Model.hpp"
#include "safepll.h"


Model::Model(unsigned int states, unsigned int rate_cats)
{
  set_simple_model(states, rate_cats);
}

void Model::set_simple_model(unsigned int states, unsigned int rate_cats)
{
  frequencies = std::vector<double>(states, 1.0/double(states));
  substitution_parameters = std::vector<double> ( (states * (states - 1)) / 2, 1.0);
  rate_categories.resize(rate_cats);
  pll_compute_gamma_cats(1, rate_cats, &rate_categories[0]);
  params_indices = std::vector<unsigned int> (rate_cats, 0);
}


