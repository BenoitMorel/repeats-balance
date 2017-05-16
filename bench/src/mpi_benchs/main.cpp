
#include "mpi_benchs.hpp"
#include <iostream>

int main(int argc, char *params[])
{
  std::cout << "Main " << std::endl;
  if (argc < 2) {
    std::cerr << "Error : missing parameter" << std::endl;
    return 1;
  }
  if (!strcmp(params[1], "full_tree_traversal")) {
    full_tree_traversal(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }



  return 0;
}


