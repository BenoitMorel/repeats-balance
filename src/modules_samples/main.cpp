
#include "modules_samples.h"

#include <iostream>

int main(int argc, char *params[])
{
  
  if (argc < 2) {
    std::cerr << "Error : missing parameter" << std::endl;
    return 1;
  }
  if (!strcmp(params[1], "terraces_experiments")) {
    terraces_experiments(argc - 2, params + 2);
  } else {
    std::cerr << "Error unkown parameter " << params[1] << std::endl;
  }



  return 0;
}
