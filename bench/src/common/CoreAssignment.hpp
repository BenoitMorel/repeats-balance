#ifndef CORE_ASSIGNMENT_HH
#define CORE_ASSIGNMENT_HH

#include <vector>
#include "PartitionIntervals.hpp"

class CoreAssignment {


public:
  CoreAssignment(unsigned int core_id) : core_id(core_id), weight(0) {}

  void add_assignment(unsigned int msa_id, unsigned int start, unsigned int size, double weight) {
    this->weight += weight;
    PartitionIntervals interval(msa_id);
    interval.add_interval(start, size);
    assignment.push_back(interval);
  }

  unsigned int get_core_id() const {return core_id;}
  
  double get_core_weight() const {return weight;}

  friend std::ostream& operator<< (std::ostream& stream, const CoreAssignment& core_as) {
    stream << "{";
    for (unsigned int i = 0; i < core_as.assignment.size(); ++i) {
      stream << core_as.assignment[i] << ", ";
    }
    stream << "}";
    return stream;
  }

  const std::vector<PartitionIntervals> &get_assignments() const {return assignment;}
private:
  unsigned int core_id;
  std::vector<PartitionIntervals> assignment;
  double weight;
};


#endif
