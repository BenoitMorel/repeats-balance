#ifndef _RB_PARTITION_
#define _RB_PARTITION_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

class Partition; 
class InputSequences;

typedef std::vector<Partition> Partitions;
typedef std::vector<const Partition *> PartitionsPointers;



class Partition {
  public:
    Partition() : _sequences(0) {

    }

    void init(const InputSequences *sequences, unsigned int index, unsigned int start, unsigned int size) {
      _sequences = sequences;
      _start = start;
      _site_costs.resize(size);
      std::fill(_site_costs.begin(), _site_costs.end(), 0.0);
      _index = index;
    }

    unsigned int index() const {
      return _index;
    }
    
    unsigned int start() const {
      return _start;
    }

    unsigned int size() const {
      return _site_costs.size();
    }

    const InputSequences *sequences() const {
      return _sequences;
    }

    std::vector<double> &site_costs() {
      return _site_costs;
    }
     
    // for tests
    void fill_costs_randomly() {
      for (unsigned int i = 0; i < _site_costs.size(); ++i) {
        _site_costs[i] = double(rand()) / double(RAND_MAX);
      }
    } 

    void normalize_costs(double ratio) {
      _accumulated.resize(_site_costs.size());
      _site_costs[0] = 1 - _site_costs[0] / ratio;
      _accumulated[0] = _site_costs[0];
      for (unsigned int i = 1; i < size(); ++i) {
        _site_costs[i] = 1 - _site_costs[i] / ratio;
        _accumulated[i] = _accumulated[i - 1] + _site_costs[i];
      } 
    }

    // only makes sense after normalize_costs was called
    double total_weight() const {
      return _accumulated[_accumulated.size() - 1];
    }

    unsigned int find_upper_bound(double value) const {
      return  std::upper_bound(_accumulated.begin(), _accumulated.end(), value) - _accumulated.begin();
    }

    double get_accumulated_weight(unsigned int site) const {
      return _accumulated[site];
    }

    friend std::ostream& operator<< (std::ostream &out, const Partition &part) {
      out << "Partition starting from " << part.start() << std::endl;
      out << "size : " << part.size() << std::endl;
      out << "average site repeats ratio : ";
      for (unsigned int i = 0; i < part.size(); ++i) {
        out << part._site_costs[i] << " ";
      }
      out << std::endl;
      return out;
    }

 
	static void get_sorted_partitions(const Partitions &partitions, PartitionsPointers &o_partitions_ptr) {
    o_partitions_ptr.resize(partitions.size());
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      o_partitions_ptr[i] = &partitions[i];
    }
    std::sort(o_partitions_ptr.begin(), o_partitions_ptr.end(), compare_ptr_partitions);
  }
	
  static void get_sorted_partitions_weighted(const Partitions &partitions, PartitionsPointers &o_partitions_ptr) {
    o_partitions_ptr.resize(partitions.size());
    for (unsigned int i = 0; i < partitions.size(); ++i) {
      o_partitions_ptr[i] = &partitions[i];
    }
    std::sort(o_partitions_ptr.begin(), o_partitions_ptr.end(), compare_ptr_partitions_weight);
  }

  private:
	  static bool compare_ptr_partitions(const Partition* a, const Partition* b) { return (a->size() < b->size()); }
	  static bool compare_ptr_partitions_weight(const Partition* a, const Partition* b) { return (a->total_weight() < b->total_weight()); }
    unsigned int _index;
    unsigned int _start;
    const InputSequences *_sequences;
    std::vector<double> _site_costs;
    std::vector<double> _accumulated;
};



#endif
