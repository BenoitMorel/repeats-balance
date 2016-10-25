#ifndef _SB_TREE_H_
#define _SB_TREE_H_

#include <vector>
#include <algorithm>
#include "repeatsbalance.h"
#include "Node.hpp"
class Tree {
  public:
    Tree() : _buffer(100000000) {
      std::fill(_buffer.begin(), _buffer.end(), 0); 
    }

    /* 
     * generate a random tree with leaves_number leaves
     */
    void set_random(unsigned int leaves_number) {
      reset();
      _root = Node::generate_random_tree(leaves_number, _nodes_pool);
      _root->get_postorder(_post_order_nodes);
    }

    void update_SRcount(Partition &partition) {
      update_SRcount(*partition.sequences(), partition.start(), partition.site_costs());
    }

    /*
     *  Update the site repeats count vector. Each element of the array corresponds to a site. The method
     *  adds one to this element for each new repeat in this element according to this tree.
     */
    void update_SRcount(const InputSequences &sequences, unsigned int offset,  std::vector<double> &o_SRCount) {
      _cleanbuffer.resize(sequences.width());
      for (unsigned int i = 0; i < _post_order_nodes.size(); ++i) {
        _post_order_nodes[i]->fill_identifier(sequences, offset, _buffer, _cleanbuffer,  o_SRCount);
      }
    }

    /**
     *  log
     */
    friend std::ostream& operator<< (std::ostream &out, const Tree &tree) {
      out << *tree._root;
      return out;
    }
  private:
    
    void reset() {
      _root = 0;
      _nodes_pool.clear();
      _post_order_nodes.clear();
    }

    Node *_root;
    std::vector<Node> _nodes_pool;
    std::vector<Node *> _post_order_nodes;
    std::vector<int> _buffer;
    std::vector<int> _cleanbuffer;
};

#endif

