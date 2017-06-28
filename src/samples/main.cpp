
#include "samples.hpp"
#include <iostream>
#include <cstring>

int main(int argc, char *params[])
{
  
  if (argc < 2) {
    std::cerr << "Error : missing parameter" << std::endl;
    return 1;
  }
  if (!strcmp(params[1], "print_random_trees")) {
    print_random_trees(argc - 2, params + 2);
  } else if (!strcmp(params[1], "print_all_trees")) {
    print_all_trees(argc - 2, params + 2);
  } else if (!strcmp(params[1], "random_trees_likelihoods")) {
    random_trees_likelihoods(argc - 2, params + 2);
  } else if (!strcmp(params[1], "tbnni_move")) {
    tbnni_move(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }



  return 0;
}

