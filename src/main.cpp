#include <iostream>
#include <vector>
#include "Node.hpp"

int check(std::vector<Node> &nodes, int seq_size, int *expectedSRNumbers, const char *test_name) {
  std::vector<int> srcounts(seq_size);
  for (unsigned int node = 0; node < nodes.size(); ++node) {
    nodes[node].id = node;
    nodes[node].fill_identifier(seq_size, srcounts);
  } 
  std::cout << "srcounts : " ;
  for (unsigned int i = 0 ; i < srcounts.size(); ++i) {
    if (srcounts[i] != expectedSRNumbers[i]) {
      std::cout << "ERROR IN " << test_name << std::endl;
      std::cout << "Expected: "; 
      for (int j = 0; j < seq_size; ++j) {
        std::cout << expectedSRNumbers[j] << ' ';
      }
      std::cout << std::endl << "Found:    ";
      for (int j = 0; j < seq_size; ++j) {
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
  const int seq_number = 4;
  const int seq_size = 4;
  const int nodes_number = 7;
  const char seq[seq_number][seq_size] = 
                    {{'A','A','A','C'},     //   0---2------6 
                     {'C','C','C','C'},     //   1---'      '
                     {'G','G','G','G'},     //   3---5------'
                     {'T','T','A','A'}};    //   4---'

  std::vector<Node> nodes(nodes_number);
  nodes[6].set_children(&nodes[2], &nodes[5]);
  nodes[2].set_children(&nodes[0], &nodes[1]);
  nodes[5].set_children(&nodes[3], &nodes[4]);
  nodes[0].set_sequence(seq[0]);
  nodes[1].set_sequence(seq[1]);
  nodes[3].set_sequence(seq[2]);
  nodes[4].set_sequence(seq[3]);
  int expectedSRNumbers[seq_size] = {0, 3, 1, 1};
  return check(nodes, seq_size, expectedSRNumbers, "test1");
}

int test2() {
  const int seq_number = 5;
  const int seq_size = 6;
  const int nodes_number = 9;
  const char seq[seq_number][seq_size] = 
                    {{'C','A','A','A','C','C'},  // 0---2-------4-----8   
                     {'A','A','G','A','A','A'},  // 1---'       '     '
                     {'A','C','G','T','T','T'},  // 3-----------'     '
                     {'T','A','A','A','T','T'},  // 5---7--------------  
                     {'G','G','T','C','G','G'}}; // 6---'  
  std::vector<Node> nodes(nodes_number);
  nodes[8].set_children(&nodes[4], &nodes[7]);
  nodes[4].set_children(&nodes[2], &nodes[3]);
  nodes[7].set_children(&nodes[5], &nodes[6]);
  nodes[2].set_children(&nodes[0], &nodes[1]);
  nodes[0].set_sequence(seq[0]);
  nodes[1].set_sequence(seq[1]);
  nodes[3].set_sequence(seq[2]);
  nodes[5].set_sequence(seq[3]);
  nodes[6].set_sequence(seq[4]);
  int expectedSRNumbers[seq_size] = {0, 0, 0, 1, 2, 4};
  return check(nodes, seq_size, expectedSRNumbers, "test2");
}



int main()
{
  test1();
  test2();
  return 0;
}
