#ifndef _RB_PARTITION_
#define _RB_PARTITION_

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
  private:
    unsigned int _start;
    std::vector<double> _site_costs;
};

typedef std::vector<Partition> Partitions;

#endif
