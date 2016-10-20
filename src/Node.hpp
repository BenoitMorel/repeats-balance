#ifndef _RB_NODE_
#define _RB_NODE_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "repeatsbalance.h"


class Node {
  public:
    Node() : _left(0), _right(0), _seq_index(0), _max_identifier(0) {      
    }

    void set_children(Node *left, Node *right) {
      _left = left;
      _right = right;
    }

    void set_sequence(int seq_index) {
      _seq_index = seq_index;
    }
    
    /*
     * Fill the site repeats identifiers in this from the input sequences and update the global o_srcounts
     * Must be called in postprocess order
     */
    void fill_identifier(const InputSequences &sequences, std::vector<int> &o_srcounts) {
      unsigned int sites_number = sequences._seq_size;
      _identifiers.resize(sites_number);
      if (!_left) { // leaf node
        // TODO dirty, must be optimized
        std::vector<int> map(256);
        std::fill(map.begin(), map.end(), 0);
        for (int site = 0; site < sites_number; site++) {
          char c = sequences._sequences[_seq_index][site];
          if (!map[c]) {
            map[c] = ++_max_identifier;
          }
          _identifiers[site] = map[c];
        }

      } else { // internal node
        // TODO allocate map only once
        std::vector<int> map(_left->_max_identifier * _right->_max_identifier); 
        std::fill(map.begin(), map.end(), 0);
        for (int site = 0; site < sites_number; ++site) {
          int index_map = _left->_identifiers[site] - 1 + _left->_max_identifier * (_right->_identifiers[site] - 1);
          if (!map[index_map]) { // new identifier
            map[index_map] = ++_max_identifier;
          } else {
            o_srcounts[site]++;
          }
          _identifiers[site] = map[index_map];
        }
      }
    }

    // not really random with respect to left and right. but it doesnt matter here
    static Node *generate_random_tree(unsigned int leafs_number, unsigned int seed, std::vector<Node> &nodes) {
      srand(seed);
      nodes = std::vector<Node>(2 * leafs_number - 1);
      std::vector<int> leaves(leafs_number);
      for (unsigned int i = 0; i < leafs_number; ++i) {
        leaves[i] = i;
      }
      std::random_shuffle(leaves.begin(), leaves.end());
      for (unsigned int i = 0; i < leafs_number; ++i) {
        nodes[i]._seq_index = leaves[i];
      }
      
      //build the initial tree
      //nodes[leafs_number] will always stay the root 
      nodes[leafs_number]._left = &nodes[0];
      nodes[leafs_number]._right = &nodes[1];

      //stepwise addition
      for (unsigned int i = 2; i < leafs_number; ++i) { 
        int parent = (rand() % (i - 1)) + leafs_number;
        Node *curr_inner = &nodes[i - 1 + leafs_number];
        curr_inner->_left = &nodes[i];
        if (rand() % 2) {
          curr_inner->_right = nodes[parent]._left;
          nodes[parent]._left = curr_inner;
        } else { 
          curr_inner->_right = nodes[parent]._right;
          nodes[parent]._right = curr_inner;
        }
      } 
      return &nodes[leafs_number];
    }
    
    friend std::ostream& operator<< (std::ostream &out, const Node &node) {
      if (node._left) {
        out << "(" << *node._left << "," << *node._right << ")"; 
      } else {
        out << node._seq_index;
      }
      return out;
    }
    
    void get_postorder(std::vector<Node *> &postorder) {
      if (_left) {
        _left->get_postorder(postorder);
        _right->get_postorder(postorder);
      }
      postorder.push_back(this);
    }

  private:
    Node *_left;
    Node *_right;
    unsigned int _seq_index;
    unsigned int _max_identifier;
    std::vector<int> _identifiers; // 1 to _max_identifier
};

#endif
