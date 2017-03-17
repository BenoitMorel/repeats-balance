#include "common.h"


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
#define RATE_CATS 4

  void my_update_partials(pll_partition_t * partition,
                          const pll_operation_t * operations,
                          unsigned int ops_count,
                          unsigned int update_repeats)
 {
    pll_update_partials(partition, operations, ops_count);
 }

  static void fatal(const char * format, ...) __attribute__ ((noreturn));

  typedef struct
  {
    int clv_valid;
  } node_info_t;

  /* a callback function for performing a full traversal */
  static int cb_full_traversal(pll_utree_t * node)
  {
    node->data = (void *)~0;
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



PLLHelper::PLLHelper(const char * newick,
    const char * seq,
    unsigned int attribute,
    unsigned int states)
{
  do_opdate_operations = true;
  unsigned int i;
  unsigned int sequence_count;
  bool isdna = (states == 4);
  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  tree = pll_utree_parse_newick(newick, &tip_nodes_count);
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
      isdna ? pll_map_nt : pll_map_aa,
      tip_nodes_count,
      &(msa->length));


  partition = pll_partition_create(tip_nodes_count,
      inner_nodes_count,
      states,
      (unsigned int)(msa->length),
      1,
      branch_count,
      RATE_CATS,
      inner_nodes_count,
      attribute);
  if (!partition)
  {
    fatal("failed to create partition");
  }
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[4] = {0};
  pll_compute_gamma_cats(1, 4, rate_cats);
  pll_set_frequencies(partition, 0, isdna ? frequencies : pll_aa_freqs_dayhoff);
  pll_set_subst_params(partition, 0, isdna ? subst_params : pll_aa_rates_dayhoff);
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
    pll_set_tip_states(partition, tip_clv_index, 
        isdna ? pll_map_nt : pll_map_aa,
        msa->sequence[i]);
  }
  pll_msa_destroy(msa);
  hdestroy();
  free(data);
  free(tipnodes);
  travbuffer = (pll_utree_t **)malloc(nodes_count * sizeof(pll_utree_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  memset(params_indices, 0, 4 * sizeof(unsigned int));
  
  operations = (pll_operation_t *)malloc(inner_nodes_count *
      sizeof(pll_operation_t));

  sumtable = (double *)pll_aligned_alloc(
      partition->sites * partition->rate_cats * partition->states_padded *
      sizeof(double), partition->alignment);
}

double rand_double()
{
  return ((double) rand() / (RAND_MAX)); 
}

void PLLHelper::generate_random_model()
{
  double frequencies[4] = { rand_double(), rand_double(), rand_double(), rand_double() };
  double subst_params[6] = {rand_double(), rand_double(), rand_double(), rand_double(), rand_double(), rand_double()};
  double rate_cats[4] = {0};
  pll_compute_gamma_cats(1, 4, rate_cats);
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

}

void PLLHelper::set_srlookup_size(unsigned int size)
{
#ifdef HAS_REPEATS
  pll_repeats_t *repeats = partition->repeats;
  if (!repeats) {
    return;
  }
  repeats->lookup_buffer_size = size;
  free(repeats->lookup_buffer);
  repeats->lookup_buffer = (unsigned int *)
     calloc(repeats->lookup_buffer_size, sizeof(unsigned int));
#endif
}
  
unsigned int PLLHelper::compute_attribute(unsigned int use_repeats, 
                                          unsigned int additional_attr,
                                          const char *arch)
{
  unsigned int attribute = 0;
  if (!strcmp(arch, "cpu")) {
    attribute = PLL_ATTRIB_ARCH_CPU;
  } else if (!strcmp(arch, "sse")) {
    attribute = PLL_ATTRIB_ARCH_SSE;
  } else if (!strcmp(arch, "avx")) {
    attribute = PLL_ATTRIB_ARCH_AVX;
  } else if (!strcmp(arch, "avx2")) {
    attribute = PLL_ATTRIB_ARCH_AVX2;
  } else {
    std::cerr << "Error : unknown architecture " << arch << std::endl;
    return INVALID_ATTRIBUTE;
  }
  if (use_repeats == 1) {
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  } else if(use_repeats == 0) {
    attribute |= PLL_ATTRIB_PATTERN_TIP;
  }
  attribute |= additional_attr;
  return attribute;
}

void PLLHelper::update_all_partials(bool update_repeats)
{
  unsigned int plop;
  update_partials(cb_full_traversal, update_repeats, plop);
}

void PLLHelper::update_partials(int (*cbtrav)(pll_utree_t *), 
    unsigned int update_repeats, 
    unsigned int &traversal_size)
{
  if (do_opdate_operations)
    update_operations(cbtrav, traversal_size);  
  update_partials(update_repeats);
}

void PLLHelper::update_all_repeats()
{
#ifdef HAS_REPEATS
  for (unsigned int i = 0; i < ops_count; ++i) {
    pll_update_repeats(partition, &operations[i]);
  }
#endif
}



void PLLHelper::update_operations(int (*cbtrav)(pll_utree_t *), unsigned int &traversal_size)
{
  if (!pll_utree_traverse(tree,
        cbtrav,
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
  pll_update_prob_matrices(partition,
            params_indices,
            matrix_indices,
            branch_lengths,
            matrix_count);
}


void PLLHelper::update_partials(unsigned int update_repeats)
{
  my_update_partials(partition, operations, ops_count, update_repeats);
}


void PLLHelper::update_sumtables() 
{
  pll_update_sumtable(partition,
                      tree->clv_index,
                      tree->back->clv_index,
                      PLL_SCALE_BUFFER_NONE,
                      PLL_SCALE_BUFFER_NONE,
                      params_indices, 
                      sumtable);
}
  
  
void PLLHelper::get_derivative(double *d_f, double *dd_f) 
{
    pll_update_sumtable(partition,
                      tree->clv_index,
                      tree->back->clv_index,
                      PLL_SCALE_BUFFER_NONE,
                      PLL_SCALE_BUFFER_NONE,
                      params_indices, 
                      sumtable);
  for (unsigned int i = 0; i < partition->sites; ++i) {
    std::cerr << sumtable[i] << " ";
  }
  std::cerr << std::endl;
  pll_compute_likelihood_derivatives(partition,
                                     PLL_SCALE_BUFFER_NONE,
                                     PLL_SCALE_BUFFER_NONE,
                                     42.0,
                                     params_indices,
                                     sumtable,
                                     d_f, dd_f);

}

double PLLHelper::get_likelihood() 
{
  return pll_compute_edge_loglikelihood(partition,
      tree->clv_index,
      tree->scaler_index,
      tree->back->clv_index,
      tree->back->scaler_index,
      tree->pmatrix_index,
      params_indices,
      NULL);
}
  
void PLLHelper::init_tree_stats()
{
  children_number = std::vector<unsigned int>(nodes_count);
  srclasses_number = std::vector<unsigned int>(nodes_count);
  depths = std::vector<unsigned int>(nodes_count);
  init_tree_stats_rec(tree, 0);
  init_tree_stats_rec(tree->back, 0);
}

void PLLHelper::init_tree_stats_rec(pll_utree_t *root, unsigned int depth)
{
  unsigned int index = root->clv_index;
  depths[index] = depth;
#ifdef HAS_REPEATS
  if (partition->repeats) {
    srclasses_number[index] = partition->repeats->pernode_max_id[index]; 
    srclasses_number[index] = srclasses_number[index] ? srclasses_number[index] : partition->sites;
  }
#endif
  if (!root->next) {
    children_number[index] = 0;
    return;
  }
  init_tree_stats_rec(root->next->back, depth + 1);
  init_tree_stats_rec(root->next->next->back, depth + 1);
  children_number[index] = 2 
    + children_number[root->next->back->clv_index]
    + children_number[root->next->next->back->clv_index];
}

void fill_op(const pll_utree_t &tree, pll_operation_t &op)
{
  op.parent_clv_index = tree.clv_index;
  op.parent_scaler_index = tree.scaler_index;

  op.child1_clv_index = tree.next->back->clv_index;
  op.child1_scaler_index = tree.next->back->scaler_index;
  op.child1_matrix_index = tree.next->back->pmatrix_index;

  op.child2_clv_index = tree.next->next->back->clv_index;
  op.child2_scaler_index = tree.next->next->back->scaler_index;
  op.child2_matrix_index = tree.next->next->back->pmatrix_index;
}

bool coin() {
  return rand() % 2;
}

void PLLHelper::update_partial(pll_operation_t &operation, 
      unsigned int iterations, 
      bool update_repeats) 
{
  std::vector<pll_operation_t> ops(iterations);
  std::fill(ops.begin(), ops.end(), operation);
  my_update_partials(partition, &ops[0], iterations, update_repeats);
}
void PLLHelper::update_partials(std::vector<pll_operation_t> &operations, 
      unsigned int iterations, 
      bool update_repeats)
{
  for (unsigned int i = 0; i < iterations; ++i) {
    my_update_partials(partition, &operations[0], operations.size(), update_repeats);
  }
}

  
void PLLHelper::print_op_stats(pll_operation_t &op) const
{
  std::cout 
            << srclasses_number[op.parent_clv_index] << " "
            << srclasses_number[op.child1_clv_index] << " "
            << srclasses_number[op.child2_clv_index] 
            << std::endl;
}

PLLHelper::~PLLHelper() {
  pll_partition_destroy(partition);
  //pll_utree_destroy(tree);
  pll_aligned_free(sumtable);
  free(travbuffer);
  free(branch_lengths);
  free(operations);
  free(matrix_indices);
}


Timer::Timer() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start); 
}

long Timer::get_time() {
  timespec end;
  timespec temp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return (temp.tv_sec * 1000000000 + temp.tv_nsec) / 1000000; 
}




