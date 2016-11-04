
#include <iostream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"
#include "Helper.hpp"

void experiment1() {
//  Helper::experiment1("../data/94/94.phy", "../data/94/94.part", 10, 30, "../results/experiment1/94_30cpus");
  Helper::experiment1("../data/94/94.phy", "../data/94/94.part", 10, 200, "../results/experiment1/94_200cpus");  
  Helper::experiment1("../data/59/59.phy", "../data/59/59.part", 10, 10, "../results/experiment1/59_10cpus");
  Helper::experiment1("../data/128/128.phy", "../data/128/128.part", 10, 10, "../results/experiment1/128_10cpus");
  Helper::experiment1("../data/404/404.phy", "../data/404/404.part", 10, 10, "../results/experiment1/404_10cpus");
}

void experiment2() {
  Helper::experiment2("../data/59/59.phy", "../data/59/59.part", "../data/59/seed12345_top0.newick", 
                      10, 10, "../results/experiment2/59_10cpus_seed12345_top0");
/*  Helper::experiment2("../data/128/128.phy", "../data/128/128.part", "../data/128/firsttree.newick", 
                      10, 10, "../results/experiment2/128_10cpus_first");
  Helper::experiment2("../data/128/128.phy", "../data/128/128.part", "../data/128/tree86894.newick", 
                      10, 10, "../results/experiment2/128_10cpus_tree86894");
  Helper::experiment2("../data/404/404.phy", "../data/404/404.part", "../data/404/seed42_firsttree.newick", 
                      10, 10, "../results/experiment2/404_10cpus_first");
  Helper::experiment2("../data/404/404.phy", "../data/404/404.part", "../data/404/seed42_tree_98240.newick", 
                      10, 10, "../results/experiment2/404_10cpus_98240");
*/
}

void experiment3() {
  Helper::experiment3("../data/128/split-partitions/12S_rRNA.phy", "../results/experiment3/exp3_12S_rRNA.tex");
  Helper::experiment3("../data/128/split-partitions/CREM.phy", "../results/experiment3/exp3_CREM.tex");
  Helper::experiment3("../data/128/split-partitions/BRCA1.phy", "../results/experiment3/exp3_BRCA1.tex");
}

void experiment4() {
  Helper::experiment4("../data/404/404.phy", "../data/404/404.part", "../data/404/treechar_samples/tree404.", 
                      53, 50, 10, "../results/experiment4/404_10cpus.tex");
}

void experiment5() {
  Helper::experiment5("../data/404/404.phy", "../data/404/404.part",  
                      30, 50, 10, "../results/experiment5/404_10cpus.tex");
  Helper::experiment5("../data/59/59.phy", "../data/59/59.part",  
                      30, 50, 10, "../results/experiment5/59_10cpus.tex");
  Helper::experiment5("../data/128/128.phy", "../data/128/128.part",  
                      30, 50, 10, "../results/experiment5/128_10cpus.tex");
}

int main()
{
  int seed = time(0);
  seed = 1477856748;
  std::cout << "seed : " << seed << std::endl;
  srand(seed);
  //experiment1();
  //experiment2(); 
  //experiment3(); 
  experiment4(); 
  experiment5(); 
  return 0;
}


