
#include "common.h"
#include <iostream>
int main(int argc, char *params[])
{
  
  if (argc < 2) {
    std::cerr << "Error : missing parameter" << std::endl;
    return 1;
  }
  if (!strcmp(params[1], "full_traversal")) {
    full_traversal(argc - 2, params + 2);
  } else if (!strcmp(params[1], "derivatives")) {
    derivatives(argc - 2, params + 2);
  } else if (!strcmp(params[1], "sumtables")) {
    sumtables(argc - 2, params + 2);
  } else if (!strcmp(params[1], "update_repeats")) {
    update_repeats(argc - 2, params + 2);
  } else if (!strcmp(params[1], "pernode")) {
    pernode(argc - 2, params + 2);
  } else if (!strcmp(params[1], "calibrate_test")) {
    calibrate_test(argc - 2, params + 2);
  } else if (!strcmp(params[1], "numerics")) {
    numerics(argc - 2, params + 2);
  } else if (!strcmp(params[1], "repeats_rates")) {
    repeats_rates(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }



  return 0;
}
