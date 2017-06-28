#include "samples.hpp"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "../common/repeatsbalance.hpp"


// copied from libpll
static void print_node_info(const pll_unode_t * node, int options)
{
  if (options & PLL_UTREE_SHOW_LABEL)
    printf (" %s", node->label);
  if (options & PLL_UTREE_SHOW_BRANCH_LENGTH)
    printf (" %f", node->length);
  if (options & PLL_UTREE_SHOW_CLV_INDEX)
    printf (" %d", node->clv_index);
  if (options & PLL_UTREE_SHOW_SCALER_INDEX)
    printf (" %d", node->scaler_index);
  if (options & PLL_UTREE_SHOW_PMATRIX_INDEX)
    printf (" %d", node->pmatrix_index);
  printf("\n");
}

// copied from libpll
static void utree_link(pll_unode_t * a,
    pll_unode_t * b,
    double length,
    unsigned int pmatrix_index)
{
  a->back = b;
  b->back = a;
  a->length = length;
  b->length = length;

  a->pmatrix_index = b->pmatrix_index = pmatrix_index;
}

// copied from libpll
static void utree_swap(pll_unode_t * t1, pll_unode_t * t2) 
{
  /* swaps the positions of trees t1 and t2. The two trees retain the branch
   *   lengths from their root to their respective parent nodes, and retain their
   *     pmatrix indices (i.e. no updating of pmatrices is required) */

  pll_unode_t * temp = t1->back;

  utree_link(t1, t2->back, t2->back->length, t2->back->pmatrix_index);
  utree_link(t2, temp, temp->length, temp->pmatrix_index);
}



//        A
//        |
//       / \
//      /    \
//  B--|u   v|--E
//     |     |
//     C     D
unsigned int tbnni_rollback_table[] = {0, 1, 2, 3, 4, 5, 6, 10, 14, 9, 7, 11, 12, 13, 8};
struct tbnni_rb {
  pll_unode_t *node;
  unsigned int tbnni_number;
};

void apply_tbnni_move(pll_unode_t *unode, unsigned int tbnni_number, tbnni_rb *rb) {
  unsigned int subtree_move =  tbnni_number % 3;
  unsigned int singlenode_move =  tbnni_number / 3;
  pll_unode_t *nodes[5];
  pll_unode_t *u = unode->next->back;
  pll_unode_t *v = unode->next->next->back;
  nodes[0] = unode;         // A
  nodes[1] = u->next;       // B
  nodes[2] = u->next->next; // C
  nodes[3] = v->next;       // D
  nodes[4] = v->next->next; // E
  if (singlenode_move) {
    utree_swap(nodes[0], nodes[singlenode_move]);
  }
  if (subtree_move) {
    utree_swap(nodes[2], nodes[2 + subtree_move]); // C and D
  }
  if (rb) {
    rb->node = nodes[0];
    rb->tbnni_number = tbnni_rollback_table[tbnni_number];
  }
}




void tbnni_move(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip part states newick" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  const char *newick_filename = params[i++];

  srand(time(0));
  MSA msa(phy_filename, states_number);
  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  std::cout.precision(17);
  for (unsigned int i = 0; i < 15; ++i) {
    std::cout << i << std::endl;
    Tree tree(&msa, newick_filename);
    //tree.print();
    pll_utree_t *utree = tree.get_pll_tree();
    pll_unode_t *unode = utree->nodes[0]->back->next->next->back;
    tbnni_rb rb;
    apply_tbnni_move(unode, i, &rb);
    tree.print();
    apply_tbnni_move(rb.node, rb.tbnni_number, 0);
    tree.print();
    std::cout << std::endl;
    /*
    LikelihoodEngine engine(&tree, &msa, part_filename, attribute, states_number, 4, 0);  
    engine.update_operations();
    engine.update_matrices();
    engine.update_partials();
    std::cout << engine.compute_likelihood() << std::endl;
    */
  }

}


