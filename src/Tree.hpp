#ifndef _SB_TREE_H_
#define _SB_TREE__

#include <vector>
#include <algorithm>
#include "repeatsbalance.h"
#include "Node.hpp"
class Tree {
  public:
    Tree() : _buffer(100000000) {
      std::fill(_buffer.begin(), _buffer.end(), 0); 
    }


    void set_random(unsigned int leaves_number, unsigned int seed) {
      reset();
      _root = Node::generate_random_tree(leaves_number, seed, _nodes_pool);
      _root->get_postorder(_post_order_nodes);
    }

    void reset() {
      _root = 0;
      _nodes_pool.clear();
      _post_order_nodes.clear();
    }

    void update_SRcount(const InputSequences &sequences, unsigned int offset,  std::vector<double> &o_SRCount) {
      _cleanbuffer.resize(sequences.width());
      for (unsigned int i = 0; i < _post_order_nodes.size(); ++i) {
        _post_order_nodes[i]->fill_identifier(sequences, offset, _buffer, _cleanbuffer,  o_SRCount);
      }
    }

    friend std::ostream& operator<< (std::ostream &out, const Tree &tree) {
      out << *tree._root;
      return out;
    }
  private:
    Node *_root;
    std::vector<Node> _nodes_pool;
    std::vector<Node *> _post_order_nodes;
    std::vector<int> _buffer;
    std::vector<int> _cleanbuffer;
};

#endif

