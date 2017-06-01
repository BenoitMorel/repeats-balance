
#include "tools.hpp"
#include <iostream>
#include <cstring>

int main(int argc, char *params[])
{
  
  if (argc < 2) {
    std::cerr << "Error : missing parameter" << std::endl;
    return 1;
  }
  if (!strcmp(params[1], "analyse_partitions")) {
    analyse_partitions(argc - 2, params + 2);
  } else if (!strcmp(params[1], "generate_repeats_dataset")) {
    generate_repeats_dataset(argc - 2, params + 2);
  } else if (!strcmp(params[1], "generate_partitions")) {
    generate_partitions(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }



  return 0;
}

