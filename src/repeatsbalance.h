#ifndef _REPEATBALANCE_H_
#define _REPEATBALANCE_H_

#define SRLOG_ENABLED

#ifdef SRLOG_ENABLED
  #include <iostream>
  #define SRLOG(msg) \
      std::cout << __FILE__ << "(" << __LINE__ << "): " << msg << std::endl 
#else
  #define SRLOG(msg)
#endif

#include <string>

class InputSequences {
  public:
    ~InputSequences() {
      for (unsigned int i = 0; i < _seq_number; ++i) {
        delete[] _names[i];
        delete[] _sequences[i];
      } 
    }

    void add_sequence(char *name, char *sequence) {
      _names.push_back(name);
      _sequences.push_back(sequence);
    }

    unsigned int _seq_number;
    unsigned int _seq_size;
    std::vector<char *> _names;
    std::vector<char *> _sequences;

    friend std::ostream& operator<< (std::ostream &out, const InputSequences &seq) {
      out << seq._seq_number << " " << seq._seq_size << std::endl;
      for (unsigned int i = 0; i < seq._seq_number; ++i) {
        out << seq._names[i] << " " << seq._sequences[i] << std::endl;
      }
      return out;
    }

};

class InputPartitions {
  public:
    void add_partition(const std::string &name, unsigned int offset, unsigned int size) {
      _names.push_back(name);
      _offsets.push_back(offset);
      _sizes.push_back(size);
    }

    const std::string name(unsigned int i) const {
      return _names[i];
    }

    unsigned int offset(unsigned int i) const {
      return _offsets[i];
    }

    unsigned int size(unsigned int i) const {
      return _sizes[i];
    }

    unsigned int size() const {
      return _names.size();
    }

    friend std::ostream& operator<< (std::ostream &out, const InputPartitions &part) {
      for (unsigned int i = 0; i < part._names.size(); ++i) {
        out << part._names[i] << " " << part._offsets[i] << " " << part._sizes[i] << std::endl;
      }

      return out;
    }

  private:
    std::vector<std::string> _names;
    std::vector<unsigned int> _offsets;
    std::vector<unsigned int> _sizes;
};

void parse_sequences(const char *file, InputSequences &sequences);
void parse_partitions(const char *file, InputPartitions &sequences);

#endif

