#ifndef _SB_TREE_H_
#define _SB_TREE__

#include <vector>
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

    unsigned int get_nodes_number() const {
      return _nodes_pool.size();
    }

    Node *get_node(unsigned int i) {
      return _post_order_nodes[i];
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

