#ifndef _RB_NODE_
#define _RB_NODE_

#include <iostream>
#include <vector>



class Node {
  public:
    Node() : _left(0), _right(0), _seq(0), _max_identifier(0) {      
    }

    void set_children(Node *left, Node *right) {
      _left = left;
      _right = right;
    }

    void set_sequence(const char *seq) {
      _seq = seq;
    }
    
    /*
     * Fills the site repeats identifiers in this and update the global o_srcounts
     * Must be called in postprocess order
     */
    void fill_identifier(int sitesNumber, std::vector<int> &o_srcounts) {
      _identifiers.resize(sitesNumber);
      if (_seq) { // leaf node
        // TODO dirty, must be optimized
        std::vector<int> map(256);
        std::fill(map.begin(), map.end(), 0);
        for (int site = 0; site < sitesNumber; site++) {
          if (!map[_seq[site]]) {
            map[_seq[site]] = ++_max_identifier;
          }
          _identifiers[site] = map[_seq[site]];
        }

      } else { // internal node
        // TODO allocate map only once
        std::vector<int> map(_left->_max_identifier * _right->_max_identifier); 
        std::fill(map.begin(), map.end(), 0);
        for (int site = 0; site < sitesNumber; ++site) {
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
    
    friend std::ostream& operator<< (std::ostream &out, const Node &node) {
      out << "id : " << node.id << std::endl;
      out << "identifiers : "; 
      for (unsigned int i = 0 ; i < node._identifiers.size(); ++i) {
        out << node._identifiers[i] << ' ';
      }
      out << std::endl;
      return out;
    }

    int id;
    


  private:
    Node *_left;
    Node *_right;
    const char *_seq;
    int _max_identifier;
    std::vector<int> _identifiers; // 1 to _max_identifier
};

#endif
