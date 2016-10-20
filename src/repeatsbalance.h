#ifndef _REPEATBALANCE_H_
#define _REPEATBALANCE_H_


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

void parse_sequences(const char *file, InputSequences &sequences);

#endif

