
#include "libpll_helper.h"
#include "pll.h"
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
#include <iostream>
#define STATES    4
#define RATE_CATS 4
#ifndef PLL_ATTRIB_SITES_REPEATS // old pll version
  #define PLL_ATTRIB_SITES_REPEATS 0
#endif
  static void fatal(const char * format, ...) __attribute__ ((noreturn));

  typedef struct
  {
    int clv_valid;
  } node_info_t;

  /* a callback function for performing a full traversal */
  static int cb_full_traversal(pll_utree_t * node)
  {
    return 1;
  }


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

static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}



double compute_partition(const char * newick,
    const char * seq,
    unsigned int attribute,
    unsigned int iterations)
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int sequence_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_utree_t ** travbuffer;

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(newick, &tip_nodes_count);
  if (!tree)
    fatal("Tree must be an unrooted binary tree");

  /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
  set_missing_branch_length(tree, 0.000001);

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  /*  obtain an array of pointers to tip nodes */
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count,
      sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);
  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);
  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
      sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = tipnodes[i]->clv_index;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  /* read PHYLIP alignment */
  pll_msa_t * msa = pll_phylip_parse_msa(seq, &sequence_count);
  if (!msa)
    fatal(pll_errmsg);

  /* compress site patterns */
  if (sequence_count != tip_nodes_count)
    fatal("Number of sequences does not match number of leaves in tree");

  unsigned int * weight = pll_compress_site_patterns(msa->sequence,
      pll_map_nt,
      tip_nodes_count,
      &(msa->length));


  partition = pll_partition_create(tip_nodes_count,
      inner_nodes_count,
      STATES,
      (unsigned int)(msa->length),
      1,
      branch_count,
      RATE_CATS,
      inner_nodes_count,
      attribute);

  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[4] = {0};
  pll_compute_gamma_cats(1, 4, rate_cats);
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);
  pll_set_pattern_weights(partition, weight);
  free(weight);
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;
    found = hsearch(query,FIND);
    if (!found)
      fatal("Sequence with header %s does not appear in the tree", msa->label[i]);
    unsigned int tip_clv_index = *((unsigned int *)(found->data));
    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, msa->sequence[i]);
  }
  pll_msa_destroy(msa);
  hdestroy();
  free(data);
  free(tipnodes);
  travbuffer = (pll_utree_t **)malloc(nodes_count * sizeof(pll_utree_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
      sizeof(pll_operation_t));
  unsigned int traversal_size;
  if (!pll_utree_traverse(tree,
        cb_full_traversal,
        travbuffer,
        &traversal_size))
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");
  pll_utree_create_operations(travbuffer,
      traversal_size,
      branch_lengths,
      matrix_indices,
      operations,
      &matrix_count,
      &ops_count);
  unsigned int params_indices[4] = {0,0,0,0};
  pll_update_prob_matrices(partition,
      params_indices,
      matrix_indices,
      branch_lengths,
      matrix_count);
  double logl = 0;
  pll_update_partials(partition, operations, ops_count);
  logl = pll_compute_edge_loglikelihood(partition,
      tree->clv_index,
      tree->scaler_index,
      tree->back->clv_index,
      tree->back->scaler_index,
      tree->pmatrix_index,
      params_indices,
      NULL);
  for (i = 1; i < iterations; ++i) 
  {
#if PLL_ATTRIB_SITES_REPEATS != 0
    pll_update_partials_top(partition, operations, ops_count, !(attribute & (1 << 10)));
#else
    pll_update_partials(partition, operations, ops_count);
#endif
    logl = pll_compute_edge_loglikelihood(partition,
        tree->clv_index,
        tree->scaler_index,
        tree->back->clv_index,
        tree->back->scaler_index,
        tree->pmatrix_index,
        params_indices,
        NULL);
  }
  std::cout << "ll " << logl << std::endl;
  pll_partition_destroy(partition);
  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);
  pll_utree_destroy(tree);
  return logl;
}


double compute_partitions(const char * newick,
    const char * seqdir,
    unsigned int attribute,
    unsigned int iterations,
    std::vector<double> &times)
{
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);
  char buf[500];
  DIR *mydir;
  struct dirent *myfile;
  mydir = opendir(seqdir);
  double res = 0.0;
  while((myfile = readdir(mydir)) != NULL)
  {
    const char * seq = myfile->d_name;
    if (seq[strlen(seq) - 1] == 'y') {
      sprintf(buf, "%s/%s", seqdir, myfile->d_name);
      res += compute_partition(newick, buf, attribute, iterations);
    }
  }
  closedir(mydir);
  gettimeofday(&t2, NULL);
  double elapsed = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsed += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  times.push_back(elapsed);
  return res;
}

double compute_loadbalancing(const char * newick,
    const char * lbdir,
    unsigned int attribute,
    unsigned int iterations,
    std::vector<double> &times)
{
  char buf[500];
  DIR *mydir;
  struct dirent *myfile;
  mydir = opendir(lbdir);
  double res = 0.0;
  while((myfile = readdir(mydir)) != NULL)
  {
    const char * seqdir = myfile->d_name;
    if (seqdir[0] == 'c') {
      sprintf(buf, "%s/%s", lbdir, myfile->d_name);
      res += compute_partitions(newick, buf, attribute, iterations, times);
    }
  }
  closedir(mydir);
  return res;
}


