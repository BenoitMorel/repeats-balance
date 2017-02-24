#ifndef PARTITION_BUILDER_H
#define PARTITION_BUILDER_H

#include "safepll.h"
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace CA {

class Timer {
  public:
    Timer() {
      reset();
    }

    void reset() {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    }

    long get_ns() {
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
      return (temp.tv_sec * 1000000000 + temp.tv_nsec);
    }

    long get_ms() {
      return get_ns() / 1000000;
    }
  private:
    timespec start;

};

class Calibrator {
  public:

  static pll_partition_t *create_partition(unsigned int sites,
    const std::vector<char> &states,
    unsigned int left_sites,
    unsigned int right_sites,
    unsigned int attribute,
    pll_operation_t &root_op)
  {
    std::vector<pll_operation_t> operations;
    pll_partition_t *partition = 0;
    unsigned int left_tips = log2(left_sites - 1) + 1;
    unsigned int right_tips = log2(right_sites - 1) + 1;
    unsigned int tips = left_tips + right_tips;
    unsigned int rates = 4;
    double rate_cats[4] = {0.13695378267140107,
                               0.47675185617665189,
                               0.99999999997958422,
                               2.38629436117236260};
    std::vector<double> branches(2 * (tips + 1));
    std::vector<unsigned int> matrix_indices(branches.size());
    std::vector<double> subst_params(sites, 1.0);
    unsigned int params_indices[4] = {0,0,0,0};
    for (unsigned int b = 0; b < branches.size(); ++b) {
      branches[b] = double(rand()) / double(RAND_MAX);
      matrix_indices[b] = b;
    }
    std::vector<double> freq(states.size(), 1.0/double(states.size()));
    partition = pll_partition_create(tips, tips - 1, states.size(), 
        sites, 1, branches.size(), rates, tips - 1, attribute);

    pll_set_frequencies(partition, 0, &freq[0]);
    pll_set_subst_params(partition, 0, &subst_params[0]);
    pll_set_category_rates(partition, rate_cats);


    std::vector<char> tips_chars(sites);
    char c[2];
    // left tips
    c[0] = states[0];
    c[1] = states[1];
    for (unsigned int t = 0; t < left_tips; ++t) {
      for (unsigned int s = 0; s < sites; ++s) {
        tips_chars[s] = c[(((s % left_sites) & (1 << t)) != 0)];
      }
      pll_set_tip_states(partition, t, pll_map_nt, &tips_chars[0]);
    }
    // right tips
    c[0] = states[2];
    c[1] = states[3];
    for (unsigned int t = 0; t < right_tips; ++t) {
      for (unsigned int s = 0; s < sites; ++s) {
        tips_chars[s] = c[((((s < right_sites ? s : 0)) & (1 << t)) != 0)];
      }
      pll_set_tip_states(partition, t + left_tips, pll_map_nt, &tips_chars[0]);
    }

    pll_update_prob_matrices(partition,
        params_indices,
        &matrix_indices[0],
        &branches[0],
        matrix_indices.size());


    operations.resize(tips - 1);
    // left_operations
    for (unsigned int op = 0; op < left_tips - 1; ++op) {
      operations[op].parent_clv_index = tips + op;
      operations[op].child1_clv_index = op ? tips + op - 1: 0;
      operations[op].child2_clv_index = op + 1;
      operations[op].parent_scaler_index = op;
      operations[op].child1_scaler_index = op ? op - 1 : PLL_SCALE_BUFFER_NONE;
      operations[op].child2_scaler_index = PLL_SCALE_BUFFER_NONE;
      operations[op].child1_matrix_index = op * 2;
      operations[op].child2_matrix_index = op * 2 + 1;
    }

    // right_operations
    for (unsigned int op = 0; op < right_tips - 1; ++op) {
      operations[op + left_tips - 1].parent_clv_index = tips + left_tips + op - 1;
      operations[op + left_tips - 1].child1_clv_index = op ? tips + op + left_tips - 2: left_tips;
      operations[op + left_tips - 1].child2_clv_index = op + 1 + left_tips;
      operations[op + left_tips - 1].parent_scaler_index = op + left_tips - 1;
      operations[op + left_tips - 1].child1_scaler_index = op ? op + left_tips - 2 : PLL_SCALE_BUFFER_NONE;
      operations[op + left_tips - 1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;
      operations[op + left_tips - 1].child1_matrix_index = (left_tips + op) * 2;
      operations[op + left_tips - 1].child2_matrix_index = (left_tips + op) * 2 + 1;
    }
    // root operations
    operations[tips - 2].parent_clv_index = tips + tips - 2;  
    operations[tips - 2].parent_scaler_index = tips - 2;  
    operations[tips - 2].child1_clv_index = tips + left_tips - 2;  
    operations[tips - 2].child2_clv_index = 2 * tips - 3;  
    operations[tips - 2].child1_scaler_index = left_tips - 2;  
    operations[tips - 2].child2_scaler_index = tips - 3;  
    operations[tips - 2].child1_matrix_index = tips * 2;  
    operations[tips - 2].child1_matrix_index = tips * 2 + 1;  

    pll_update_partials(partition, &operations[0], operations.size() - 1);
   
    root_op = operations[tips - 2];
    return partition;
  }

  static void print_stats(pll_partition_t *partition,
      pll_operation_t *operation)
  {
    unsigned int parent = operation->parent_clv_index;
    unsigned int left = operation->child1_clv_index;
    unsigned int right = operation->child2_clv_index;
    std::cout << "Indices " << parent << " " << left << " " << right << std::endl;
    std::cout << "sites " << partition->sites << std::endl;
    std::cout << "parent sites " << partition->repeats->pernode_max_id[parent] << std::endl;
    std::cout << "left  sites  " << partition->repeats->pernode_max_id[left] << std::endl;
    std::cout << "right sites  " << partition->repeats->pernode_max_id[right] << std::endl;
    std::cout << "tips " << partition->tips << std::endl;
  }

  static long bench_partials(pll_partition_t *partition, 
      pll_operation_t *operation,
      unsigned int iterations)
  {
    Timer t;
    for (unsigned int i = 0; i < iterations; ++i) {
      pll_update_partials(partition, operation, 1);
    }
    return t.get_ms();
  }
  
  static long bench_likelihood(pll_partition_t *partition, 
      pll_operation_t *operation,
      unsigned int iterations)
  {
    Timer t;
    unsigned int params_indices[4] = {0,0,0,0};
    for (unsigned int i = 0; i < iterations; ++i) {
      pll_compute_edge_loglikelihood(partition, 
          operation->child1_clv_index,
          operation->child1_scaler_index,
          operation->child2_clv_index,
          operation->child2_scaler_index,
          partition->tips * 2,
          params_indices,
          NULL);
    }
    return t.get_ms();
  }

  static long bench_sumtable(pll_partition_t *partition, 
      pll_operation_t *operation,
      unsigned int iterations)
  {
    double * sumtable = (double*)pll_aligned_alloc(sizeof(double) *
        partition->sites *
        partition->rate_cats *
        partition->states_padded,
        partition->alignment);
    Timer t;
    unsigned int params_indices[4] = {0,0,0,0};
    for (unsigned int i = 0; i < iterations; ++i) {
      pll_update_sumtable(partition, 
          operation->child1_clv_index,
          operation->child2_clv_index,
          operation->child1_scaler_index,
          operation->child2_scaler_index,
          params_indices,
          sumtable);
    }
    pll_aligned_free(sumtable);
    return t.get_ms();
  }


};






}

#endif
