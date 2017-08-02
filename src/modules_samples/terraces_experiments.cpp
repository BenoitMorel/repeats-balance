


#include <iostream>



extern "C" {
#include "pll_tree.h"
#include "pllmod_algorithm.h"
}
#include "../common/repeatsbalance.hpp"
/*
  printf ("Reading FASTA file: %s\n", fasta_file);

  pll_phylip_t * fd = pll_phylip_open(fasta_file, pll_map_phylip);
  if (!fd) {
    printf ("%s does not exist", fasta_file);
    return;
  }

  pll_msa_t * msa = pll_phylip_parse_interleaved(fd);
  if (!msa) {
    printf ("%s invalid", fasta_file);
    return;
  }
  pll_phylip_close(fd);

  int sites = msa->length;

  printf ("  Length of sequences: %d\n", sites);

  partition = pll_partition_create (tip_nodes_count,
                                    inner_nodes_count,
                                    STATES,
                                    (unsigned int) sites,
                                    1,
                                    branch_count,
                                    RATE_CATS,
                                    inner_nodes_count,
                                    attributes
                                    );

  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;

    found = hsearch (query, FIND);

    if (!found)
      printf ("Sequence with header %s does not appear in the tree", hdr);

    unsigned int tip_clv_index = *((unsigned int *) (found->data));

    pll_set_tip_states (partition, tip_clv_index, pll_map_nt, msa->sequence[i]);
  }

  hdestroy ();

  free (data);

  double frequencies[4] =
    { 0.25, 0.25, 0.25, 0.25 };
  pll_set_frequencies (partition, 0, frequencies);

  double subst_params[6] =
    { 1, 1, 1, 1, 1, 1 };
  pll_set_subst_params (partition, 0, subst_params);

    { 0 };
  unsigned int params_indices[4] = {0,0,0,0};
  pll_compute_gamma_cats (1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates (partition, rate_cats);

  printf ("Model paramters:\n");
  printf ("  Frequencies:  ");
  for (i = 0; i < STATES; i++)
    printf ("%.4f ", partition->frequencies[0][i]);
  printf ("\n");
  printf ("  Subst. rates: ");
  for (i = 0; i < (STATES * (STATES - 1)) / 2; i++)
    printf ("%.4f ", partition->subst_params[0][i]);
  printf ("\n");
  printf ("  Gamma rates:  ");
  for (i = 0; i < RATE_CATS; i++)
    printf ("%.4f ", partition->rates[i]);
  printf ("\n");

  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(tree,
                                                        tip_nodes_count,
                                                        1,
                                                        PLLMOD_TREE_BRLEN_LINKED);
  int params_to_optimize = 0;

  pllmod_treeinfo_init_partition(treeinfo,
                                              0,
                                              partition,
                                              params_to_optimize,
                                              0, 
                                              1.0, 
                                              params_indices,
                                              NULL  
                                              );


  double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);


  for (i = 0 ; i < 2; ++i)
  {
  }
  pllmod_treeinfo_destroy(treeinfo);
  pll_partition_destroy (partition);
//  pll_utree_destroy(treeinfo->root);

  printf ("Test OK!\n");

}
*/
void load_data(const char *phy_file,
    const char *part_file)
{
  const unsigned int states_number = 4;
  const unsigned int rates = 4;
  unsigned int attributes = PLL_ATTRIB_SITE_REPEATS | PLL_ATTRIB_ARCH_AVX;
  MSA full_msa(phy_file, states_number); 
  Tree tree(&full_msa);
  std::vector<PartitionIntervals> partitionning;
  PartitionIntervals::parse(part_file, partitionning);
  std::vector<Partition *> partitions; 
  for (unsigned int i = 0; i < partitionning.size(); ++i) 
  {
    MSA submsa(&full_msa, partitionning[i], i);
    submsa.compress();
    partitions.push_back(new Partition(&submsa,
        attributes,
        states_number,
        rates,
        0));
  }
  
  // pllmod stuff
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(tree.get_pll_root(),
                                                        tree.get_tips_number(),
                                                        partitions.size(),
                                                        PLLMOD_TREE_BRLEN_LINKED);
  for (unsigned int i = 0; i < partitionning.size(); ++i) 
  {
    unsigned int params_indices[4] = {0,0,0,0};
    const unsigned int params_to_optimize = 0;
    pllmod_treeinfo_init_partition(treeinfo,
                                  i,
                                  partitions[i]->get_partition(),
                                  params_to_optimize,
                                  0, 
                                  1.0, 
                                  params_indices,
                                  NULL  
                                  );
  }
  std::cout << "likelihood: " << pllmod_treeinfo_compute_loglh(treeinfo, 0) << std::endl;
    
}

void terraces_experiments(int argc, char *argv[])
{
  srand(42);
  load_data("../../data/404/404.phy", "../../data/404/404.part");



}

