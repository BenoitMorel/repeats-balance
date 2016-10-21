
#include <iostream>
#include <vector>
#include <time.h>
#include "repeatsbalance.h"
#include "Node.hpp"
#include "Tree.hpp"

int check(std::vector<Node> &nodes, InputSequences & sequences, int *expectedSRNumbers, const char *test_name) {
  unsigned int seq_size = sequences.width();
  std::vector<int> srcounts(seq_size);
  std::vector<int> buffer(100);
  std::fill(buffer.begin(), buffer.end(), 0);
  std::vector<int> cleanbuffer(buffer);
  for (unsigned int node = 0; node < nodes.size(); ++node) {
    nodes[node].fill_identifier(sequences, buffer, cleanbuffer, srcounts);
  } 
  std::cout << "srcounts : " ;
  for (unsigned int i = 0 ; i < srcounts.size(); ++i) {
    if (srcounts[i] != expectedSRNumbers[i]) {
      std::cout << "ERROR IN " << test_name << std::endl;
      std::cout << "Expected: "; 
      for (unsigned int j = 0; j < seq_size; ++j) {
        std::cout << expectedSRNumbers[j] << ' ';
      }
      std::cout << std::endl << "Found:    ";
      for (unsigned int j = 0; j < seq_size; ++j) {
        std::cout << srcounts[j] << ' ';
      }
      std::cout << std::endl;
      return 1;
    }
  }
  std::cout << "--- " << test_name << " ok !" << std::endl;
  return 0;
}

int test1() {
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
  return check(nodes, sequences, expectedSRNumbers, "test1");
}

int test2() {
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
  return check(nodes, sequences, expectedSRNumbers, "test2");
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
  
  if (!ok) {
    std::cout << "an error occured in test_sequences_parser" << std::endl;
    std::cout << "got : " << sequences << std::endl;
  } else {
    std::cout << "test_sequences_parser ok !" << std::endl;
  }
}

void test_partitions_parser() {
  InputPartitions partitions;
  parse_partitions("../data/128/128.part", partitions);
  bool ok = true;
  ok &= partitions.size() == 34;
  ok &= partitions.name(21) == "ATP7A";
  ok &= partitions.offset(21) == 22068;
  ok &= partitions.size(21) == 684;
  if (!ok) {
    std::cout << "an error occured in test_partitions_parser" << std::endl;
    std::cout << "got : " << partitions << std::endl;
  } else {
    std::cout << "test_partitions_parser ok !" << std::endl;
  }

}

void test_print_random_trees() {
  std::cout << "random trees :" << std::endl;
  Tree tree;
  for (int i = 0; i < 5; ++i) {
    tree.set_random(10, time(0) + i);
    std::cout << tree << std::endl;
  }
}


int main()
{
  test1();
  test2();
  test_sequences_parser();
  test_partitions_parser();
  test_print_random_trees();
  return 0;
}
