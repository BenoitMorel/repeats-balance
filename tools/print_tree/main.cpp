#include <iostream>
#include <pll.h>
#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

  static void set_missing_branch_length_recursive(pll_utree_t * tree,
      double length)
  {
    if (tree)
    {
      /* set branch length to default if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next)
    {   
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive(tree->next->back, length);
      set_missing_branch_length_recursive(tree->next->next->back, length);
    }   
  }
}
/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
}
int main(int argc, char *params[])
{
  if (argc != 3) {
    std::cout << "usage : ./print_tree path_to_newick dest.svg" << std::endl;
  }
  unsigned int tip_nodes_count;
  pll_utree_t *tree = pll_utree_parse_newick(params[1], &tip_nodes_count);  
  set_missing_branch_length(tree, 0.000001);
  pll_svg_attrib_t *plop = pll_svg_attrib_create();
  pll_utree_export_svg(tree, tip_nodes_count, plop, params[2]);
  pll_svg_attrib_destroy(plop);
	
}
