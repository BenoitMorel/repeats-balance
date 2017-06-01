#include "../../src/Helper.hpp"



/**
 *  Compute a load balancing from a given sequence, partitions, cpu numbers, and randomly
 *  generated trees
 *  Export it in an output directory (writes several subsequences files for each cpu)
 * */
void export_lbs() {
    Helper::compute_export_lb("../../data/59/59.phy", "../../data/59/59.part", 30, 10, 
                              "../../results/exports/59_10");        
    Helper::compute_export_lb("../../data/128/128.phy", "../../data/128/128.part", 30, 10, 
                              "../../results/exports/128_10");        
    Helper::compute_export_lb("../../data/404/404.phy", "../../data/404/404.part", 30, 10, 
                              "../../results/exports/404_10");        
    Helper::compute_export_lb("../../data/94/94.phy", "../../data/94/94.part", 30, 10, 
                              "../../results/exports/94_10");        
}


int main(int argc, char *argv[]) {
  export_lbs();
  return 0;
}



