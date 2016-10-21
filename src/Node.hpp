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
    void fill_identifier(const InputSequences &sequences, std::vector<int> &buffer, std::vector<int> &cleanbuffer, std::vector<int> &o_srcounts) {
      _computed = true;
      int _toclean_index = 0;
      unsigned int sites_number = sequences.width();
      std::vector<int> &map = buffer;
      //std::fill(map.begin(), map.end(), 0);
      _identifiers.resize(sites_number);
      if (!_left) { // leaf node
        for (unsigned int site = 0; site < sites_number; site++) {
          char c = sequences.sequence(_seq_index)[site];
          if (!map[c]) {
            map[c] = ++_max_identifier;
            cleanbuffer[_toclean_index++] = c;
          }
          _identifiers[site] = map[c];
        }

      } else { // internal node
        // TODO allocate map only once
        if (!_left->_computed || !_right->_computed ||_left->_max_identifier * _right->_max_identifier > map.size()) {
          // the map buffer is to small
          //SRLOG("[WARNING] Node::fill_identifier : too small buffer, ignoring this step");
          //SRLOG("[WARNING] Node::fill_identifier : rightmaxid " << _right->_max_identifier << " leftmaxid " << _left->_max_identifier);
          _computed = false;
        } else {
          for (unsigned int site = 0; site < sites_number; ++site) {
            int index_map = _left->_identifiers[site] - 1 + _left->_max_identifier * (_right->_identifiers[site] - 1);
            if (!map[index_map]) { // new identifier
              map[index_map] = ++_max_identifier;
              cleanbuffer[_toclean_index++] = index_map;
            } else {
              o_srcounts[site]++;
            }
            _identifiers[site] = map[index_map];
          }
        }
      }
      for (int i = 0; i < _toclean_index; ++i) {
        map[cleanbuffer[i]] = 0;
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
    bool _computed;
};

#endif
