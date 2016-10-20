#ifndef _SB_TREE_H_
#define _SB_TREE__

#include <vector>
#include "repeatsbalance.h"
#include "Node.hpp"
class Tree {
  public:

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

    void update_SRcount(const InputSequences &sequences, std::vector<int> &o_SRCount) {
      for (unsigned int i = 0; i < _post_order_nodes.size(); ++i) {
        _post_order_nodes[i]->fill_identifier(sequences, o_SRCount);
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
};

#endif

