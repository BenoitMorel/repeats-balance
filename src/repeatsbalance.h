#ifndef _REPEATBALANCE_H_
#define _REPEATBALANCE_H_

//#define SRLOG_ENABLED

#ifdef SRLOG_ENABLED
  #include <iostream>
  #define SRLOG(msg) \
      std::cout << __FILE__ << "(" << __LINE__ << "): " << msg << std::endl 
#else
  #define SRLOG(msg)
#endif

#include <string>
#include <cstring>
#include <map>
#include "Partition.hpp"

class Tree;

struct StrCompare : public std::binary_function<const char*, const char*, bool> {
  public:
        bool operator() (const char* str1, const char* str2) const
              { return std::strcmp(str1, str2) < 0; }
};

typedef std::map<const char*, unsigned int, StrCompare> SeqNameToInt;

class InputSequences {
  public:
    ~InputSequences() {
      for (unsigned int i = 0; i < number(); ++i) {
        delete[] _names[i];
        delete[] _sequences[i];
      } 
    }

    void add_sequence(const char *name, char *sequence) {
      if (!number()) {
        _seq_size = strlen(sequence);
      }
      _map[name] = _names.size();
      _names.push_back(name);
      _sequences.push_back(sequence);
    }

    unsigned int number() const {
      return _sequences.size();
    }

    unsigned int width() const {
      return _seq_size;
    }

    const char *sequence(unsigned int i) const {
      return _sequences[i];
    }

    const char *sequence_name(unsigned int i) const {
      return _names[i];
    }

    unsigned int get_taxa_index(const char* taxa_name) const {
      return (_map.find(taxa_name) != _map.end()) ?_map.at(taxa_name) : 0;
    }

    friend std::ostream& operator<< (std::ostream &out, const InputSequences &seq) {
      out << seq.number() << " " << seq._seq_size << std::endl;
      for (unsigned int i = 0; i < seq.number(); ++i) {
        out << seq._names[i] << " " << seq._sequences[i] << std::endl;
      }
      return out;
    }

    void write_subseq(unsigned int offset, unsigned int size, std::ostream &out) const {
      out << number() << " " << size << std::endl;
      for (unsigned int s = 0; s < number(); ++s) {
        out << _names[s] << " " << std::string(_sequences[s] + offset, size) << std::endl;
      }
    }

  private:
    unsigned int _seq_size;
    std::vector<const char *> _names;
    std::vector<const char *> _sequences;
    SeqNameToInt _map;
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

    void generate_partitions(Partitions &o_partitions, InputSequences *sequences) {
      o_partitions.resize(size()); 
      for (unsigned int i = 0; i < size(); ++i) {
        o_partitions[i].init(sequences, i, offset(i), size(i));
      }   
    }

  private:
    std::vector<std::string> _names;
    std::vector<unsigned int> _offsets;
    std::vector<unsigned int> _sizes;
  };

struct TreeScanner {
  TreeScanner(const InputSequences &iseq, Tree &itree): seq(iseq), tree(itree), current_node(0) {}
  const InputSequences &seq;
  Tree &tree;
  unsigned int current_node;
};


void parse_sequences(const char *file, InputSequences &sequences);
void parse_partitions(const char *file, InputPartitions &sequences);
void parse_tree(const char *file, const InputSequences &sequences, Tree &tree);

#endif

