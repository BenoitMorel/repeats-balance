#ifndef _LIKELIHOODENGINE_
#define _LIKELIHOODENGINE_


class Tree;
class Partition;

class LikelihoodEngine {
public:
  LikelihoodEngine(Tree &tree, Partition &partition);

  void update_operations();

  void update_matrices();

  void update_partials();

  double compute_likelihood();

  void update_sumtable();

  void compute_derivatives(double *d_f, double *dd_f); 

private:

  Tree &tree;
  Partition &partition;






};


#endif
