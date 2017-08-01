


#include <iostream>


#include <pll_tree.h>
#include <assert.h>
#include <stdarg.h>
#if(USE_HASH)
#include <search.h>
#endif
#include <time.h>

static void fatal (const char * format, ...);

/* function for printing non-binary trees */
static void print_newick_recurse(pll_unode_t * node);
static void print_newick(pll_unode_t * tree);

void terraces_experiments(int argc, char *argv[])
{
  unsigned int tree_count;

  if (argc != 3)
    fatal (" syntax: %s [trees_file] [consensus_threshold]", argv[0]);

  /* get arguments */
  char * filename = argv[1];
  double threshold = atof(argv[2]);

  /* build consensus tree */
  pll_consensus_utree_t * constree =
    pllmod_utree_consensus(filename,
                          threshold,
                          &tree_count);
  if (!constree)
    fatal("Error %d: %s\n", pll_errno, pll_errmsg);

  /* print it in NEWICK format */
  print_newick(constree->tree);

  /* clean up */
  pllmod_utree_consensus_destroy(constree);

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}


static void print_newick_recurse(pll_unode_t * node)
{
  pll_unode_t * child;

  if (pllmod_utree_is_tip(node))
  {
    printf("%s", node->label);
    return;
  }

  printf("(");
  child = node->next;
  while(child != node)
  {
    print_newick_recurse(child->back);

    if (child->next != node)
      printf(",");

    child = child->next;
  }
  printf(")");
  if (node->data)
  {
    consensus_data_t * cdata = (consensus_data_t *) node->data;
    printf("[%.3f]", cdata->support);
  }
}

static void print_newick(pll_unode_t * tree)
{
  printf("(");
  print_newick_recurse(tree->back);
  pll_unode_t * child = tree->next;
  while(child != tree)
  {
    printf(",");
    print_newick_recurse(child->back);
    child = child->next;
  }
  printf(");\n");
}

