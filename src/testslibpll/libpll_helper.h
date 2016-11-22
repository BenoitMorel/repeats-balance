#ifndef _LIBPLLHELPER_H_H_
#define _LIBPLLHELPER_H_H_


#include<vector>

double compute_loadbalancing(const char * newick,
                      const char * lbdir,
                      unsigned int attribute,
                      unsigned int iterations,
                      std::vector<double> &times);



#endif
