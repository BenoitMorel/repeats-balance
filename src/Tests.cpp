
#include <iostream>
#include <sstream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"
#include "LoadBalancing.hpp"


#define RBCHECK(value, testname, additional_error_msg) \
  if (value) {\
    std::cout << testname << " ok !" << std::endl;\
  } else {\
    std::cout << "An error occured in " << testname << std::endl;\
    std::cout << additional_error_msg << std::endl; \
  }

int check(std::vector<Node> &nodes, InputSequences & sequences, int offset, int size, int *expectedSRNumbers, const char *test_name) {
  unsigned int seq_size = size;
  std::vector<double> srcounts(seq_size);
  std::vector<int> buffer(100);
  std::fill(buffer.begin(), buffer.end(), 0);
  std::vector<int> cleanbuffer(buffer);
  for (unsigned int node = 0; node < nodes.size(); ++node) {
    nodes[node].fill_identifier(sequences, offset, buffer, cleanbuffer, srcounts);
  }
  bool ok = true; 
  for (unsigned int i = 0 ; i < srcounts.size(); ++i) {
    ok &=  (srcounts[i] == expectedSRNumbers[i]); 
  }
  RBCHECK(ok, test_name, "");
  return !ok;
}

int test_count_sr_simple_1() {
  InputSequences sequences;
  parse_sequences("../data/simple_seq/simple4-4.phy", sequences);
  const int nodes_number = 7;
  std::vector<Node> nodes(nodes_number);
  nodes[6].set_children(&nodes[2], &nodes[5]);
  nodes[2].set_children(&nodes[0], &nodes[1]);
  nodes[5].set_children(&nodes[3], &nodes[4]);
  nodes[0].set_sequence(0);
  nodes[1].set_sequence(1);
  nodes[3].set_sequence(2);
  nodes[4].set_sequence(3);
  int expectedSRNumbers[] = {0, 3, 1, 1};
  return check(nodes, sequences, 0, sequences.width(), expectedSRNumbers, "test_count_sr_simple_1");
}

int test_count_sr_simple_2() {
  InputSequences sequences;
  parse_sequences("../data/simple_seq/simple5-6.phy", sequences);
  const int nodes_number = 9;
  /*const char seq[seq_number][seq_size] = 
                    {{'C','A','A','A','C','C'},  // 0---2-------4-----8   
                     {'A','A','G','A','A','A'},  // 1---'       '     '
                     {'A','C','G','T','T','T'},  // 3-----------'     '
                     {'T','A','A','A','T','T'},  // 5---7--------------  
                     {'G','G','T','C','G','G'}}; // 6---'  
                     */
  std::vector<Node> nodes(nodes_number);
  nodes[8].set_children(&nodes[4], &nodes[7]);
  nodes[4].set_children(&nodes[2], &nodes[3]);
  nodes[7].set_children(&nodes[5], &nodes[6]);
  nodes[2].set_children(&nodes[0], &nodes[1]);
  nodes[0].set_sequence(0);
  nodes[1].set_sequence(1);
  nodes[3].set_sequence(2);
  nodes[5].set_sequence(3);
  nodes[6].set_sequence(4);
  int expectedSRNumbers[] = {0, 0, 0, 1, 2, 4};
  return check(nodes, sequences, 0, sequences.width(), expectedSRNumbers, "test_count_sr_simple_2");
}

// same sequence as before but with prefix and suffix wich will be ignored
int test_count_sr_simple_2_offset() {
  InputSequences sequences;
  parse_sequences("../data/simple_seq/simple5-6_extended.phy", sequences);
  const int nodes_number = 9;
  std::vector<Node> nodes(nodes_number);
  nodes[8].set_children(&nodes[4], &nodes[7]);
  nodes[4].set_children(&nodes[2], &nodes[3]);
  nodes[7].set_children(&nodes[5], &nodes[6]);
  nodes[2].set_children(&nodes[0], &nodes[1]);
  nodes[0].set_sequence(0);
  nodes[1].set_sequence(1);
  nodes[3].set_sequence(2);
  nodes[5].set_sequence(3);
  nodes[6].set_sequence(4);
  int expectedSRNumbers[] = {0, 0, 0, 1, 2, 4};
  return check(nodes, sequences, 6, 6, expectedSRNumbers, "test_count_sr_simple_2_offset");
}

void test_sequences_parser () {
  InputSequences sequences;
  parse_sequences("../data/minimal-6/minimal-6.phy", sequences);  
  bool ok = true;
  ok &= (6 == sequences.number());
  ok &= (60 == sequences.width());
  ok &= (std::string("Chicken") == sequences.sequence_name(1));
  ok &= (std::string("Whale") == sequences.sequence_name(5));
  ok &= (std::string("ATGGCATATCCATTCCAACTAGGTTTCCAAGATGCAGCATCACCCATCATAGAAGAGCTC") == sequences.sequence(5));
  RBCHECK(ok, "test_sequences_parser", ""); 
}

void test_partitions_parser() {
  InputPartitions partitions;
  parse_partitions("../data/128/128.part", partitions);
  bool ok = true;
  ok &= partitions.size() == 34;
  ok &= partitions.name(21) == "ATP7A";
  ok &= partitions.offset(21) == 22068;
  ok &= partitions.size(21) == 684;
  RBCHECK(ok, "test_partitions_parser", ""); 
}


void test_partitions_sort() {
  Partitions partitions(10000);
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i].init(0, i, 0, (rand() % 100000) + 1); 
  }
  PartitionsPointers partitions_ptr;
  Partition::get_sorted_partitions(partitions, partitions_ptr);
  bool ok = true;
  for (unsigned int i = 1; i < partitions.size(); ++i) {
    ok &= (partitions_ptr[i - 1]->size() <= partitions_ptr[i]->size()); 
  }
  RBCHECK(ok, "test_partitions_sort", "")
}

void test_print_random_trees() {
  std::cout << "random trees :" << std::endl;
  Tree tree;
  for (int i = 0; i < 5; ++i) {
    tree.set_random(10);
    std::cout << tree << std::endl;
  }
}

void test_naive_loadbalancing (unsigned int partitions_number, unsigned int cpus_number) {
  Partitions partitions(partitions_number);
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i].init(0, i, 0, (rand() % 100000) + 1); 
  }
  LoadBalancing lb(partitions, cpus_number);
  lb.compute_naive();
  SRLOG(lb.max_partitions_difference());
  RBCHECK(lb.is_consistent(), "test_naive_loadbalancing", "");
}

void test_kassian_loadbalancing (unsigned int partitions_number, unsigned int cpus_number) {
  Partitions partitions(partitions_number);
  for (unsigned int i = 0; i < partitions.size(); ++i) {
    partitions[i].init(0, i, 0, (rand() % 100000) + 1); 
  }
  LoadBalancing lb(partitions, cpus_number);
  lb.compute_kassian();
  bool ok = true;
  ok &= lb.is_consistent();
  ok &= lb.is_sites_balanced();
  ok &= lb.max_partitions_difference() <= 1;

  std::ostringstream os;
  os << "test_kassian_loadbalancing(" << partitions_number << "," << cpus_number << ")";
  RBCHECK(ok, os.str(), "");
}

void test_kassian_128() {
  InputPartitions ipartitions;
  parse_partitions("../data/128/128.part", ipartitions);
  Partitions partitions;
  ipartitions.generate_partitions(partitions, 0);
  bool ok = true;
  LoadBalancing lb(partitions, 4);
  lb.compute_kassian();
  ok &= lb.is_consistent();
  ok &= lb.is_sites_balanced();
  ok &= lb.max_partitions_difference() <= 1;
  RBCHECK(ok, "test_kassian_128", "");
}



int main()
{
  test_count_sr_simple_1();
  test_count_sr_simple_2();
  test_count_sr_simple_2_offset();
  test_sequences_parser();
  test_partitions_parser();
  test_partitions_sort();
  test_naive_loadbalancing(1000, 100);
  test_kassian_loadbalancing(5, 3);
  test_kassian_loadbalancing(10, 10);
  test_kassian_loadbalancing(150, 10);
  test_kassian_loadbalancing(100, 100);
  test_kassian_loadbalancing(1500, 100);
  test_kassian_128();
  return 0;
}
