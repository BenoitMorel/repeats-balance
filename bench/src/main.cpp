
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
  } else if (!strcmp(params[1], "partitioned_full_traversal")) {
    partitioned_full_traversal(argc - 2, params + 2);
  } else if (!strcmp(params[1], "pernode")) {
    pernode(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }








  return 0;
}
