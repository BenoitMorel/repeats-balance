#ifndef _RB_PARTITION_
#define _RB_PARTITION_

#include <vector>
#include <algorithm>

class Partition; 

typedef std::vector<Partition> Partitions;
typedef std::vector<const Partition *> PartitionsPointers;

class Partition {
  public:
    void init(unsigned int start, unsigned int size) {
      _start = start;
      _site_costs.resize(size);
    }
    
    unsigned int start() const {
      return _start;
    }

    unsigned int size() const {
      return _site_costs.size();
    }

    std::vector<double> &site_costs() {
      return _site_costs;
    }

    void normalize_costs(double ratio) {
      for (unsigned int i = 0; i < size(); ++i) {
        _site_costs[i] /= ratio;
      } 
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

  private:
	  static bool compare_ptr_partitions(const Partition* a, const Partition* b) { return (a->size() < b->size()); }
    unsigned int _start;
    std::vector<double> _site_costs;
};



#endif
