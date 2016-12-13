#include "../common/common.h"
#include <iostream>
#include <time.h>

static int cb_bayes(pll_utree_t * node);
int opt_ratio;
static int empty_traversal = 1;

int main(int argc, char *params[])
{
  if (argc != 9) {
    std::cout << "Expected syntax : ./rootpath_bench newick phy "
      "runs iterations use_repeats update_repeats ratio seed" << std::endl;
    return 1;
  }
  unsigned int i = 0;
  const char *tree = params[++i];
  const char *seq = params[++i];
  unsigned int runs = atoi(params[++i]);
  unsigned int iterations = atoi(params[++i]);
  bool use_repeats = atoi(params[++i]);
  bool update_repeats = atoi(params[++i]) && use_repeats;
  opt_ratio = atoi(params[++i]);
  unsigned int seed = atoi(params[++i]);
  unsigned int attribute = 0;
  srand(seed);
  if (use_repeats) {
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  } else {
    attribute |= PLL_ATTRIB_PATTERN_TIP;	
  }
  PLLHelper d(tree, seq, attribute);
  d.update_all_partials();


  for (unsigned int i = 0; i < runs; ++i) {
    unsigned int size = 100;
    empty_traversal = 1;
    d.update_operations(cb_bayes, size);
    Timer t;
    for (unsigned int j = 0; j < iterations; ++j) {
      d.update_partials(update_repeats);
    }
    unsigned int elapsed = t.get_time();
    std::cout  <<  size << "," << elapsed << std::endl;
  }
  
  return 1;
}

static int cb_bayes(pll_utree_t * node)
{
  int r;
  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next) return 0;

  if (empty_traversal == 1)
  {
    /* disable one of the directions */
    r = rand() % 3;
    switch (r)
    {
      case 0:
        node->back->data = (void *)0;
        break;
      case 1:
        node->next->back->data = (void *)0;
        break;
      case 2:
        node->next->next->back->data = (void *)0;
        break;
      default:
        assert(0);

    }
    empty_traversal = 0;
    return 1;
  }

  if (!node->data)
  {
    node->data = (void *)~0;               /* enable it for the next time we call partial traversal */
    return 0;
  }

  /* TODO: draw randomly */

  /* disable one of the directions */
  r = rand() % 3;
  switch (r)
  {
    case 0:
      node->back->data = (void *)0;
      break;
    case 1:
      node->next->back->data = (void *)0;
      break;
    case 2:
      node->next->next->back->data = (void *)0;
      break;
    default:
      assert(0);
  }
  
//  return 1;
  if ((rand() % 101) <= opt_ratio)
  {
    return 1;
  }

  return 0;
  //return (rand() % 2);
  //return 1;
}

