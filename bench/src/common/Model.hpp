#ifndef _MODEL_
#define _MODEL_

#include <vector>

struct Model {

  Model(unsigned int states, unsigned int rate_cats);
  void set_simple_model(unsigned int states, unsigned int rate_cats);

  std::vector<double> frequencies;
  std::vector<double> substitution_parameters;
  std::vector<double> rate_categories;
  std::vector<unsigned int> params_indices;


};

#endif
