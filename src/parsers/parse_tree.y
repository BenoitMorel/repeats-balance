%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "../repeatsbalance.h"
#include "../Tree.hpp"
#include "../Node.hpp"
using namespace std;


// stuff from flex that bison needs to know about:
extern  int rb_treelex();
extern  FILE *rb_treein;

void rb_treeerror(TreeScanner &scan, const char *s);
%}


%union {
  char *sval;
  Node *node;
}

%parse-param {TreeScanner &scan}

%token <sval> STRING
%token COMMA
%token RIGHT_PAR
%token LEFT_PAR
%token SEMI_COLON
%type<node> node 
%define api.prefix {rb_tree}

%%


tree: node SEMI_COLON {scan.tree.set_root($1);}
;

node: LEFT_PAR node COMMA node RIGHT_PAR {$$ = &scan.tree.get_node(scan.current_node++); $$->set_children($2, $4);}
      | LEFT_PAR node COMMA node COMMA node RIGHT_PAR {
          std::cout << "WARNING, UNROOTED TREE... adding a root" << std::endl;
          Node *node  = &scan.tree.get_node(scan.current_node++); node->set_children($2, $4);
          $$ = &scan.tree.get_node(scan.current_node++); $$->set_children(node, $6);
        }
      | STRING {$$ = &scan.tree.get_node(scan.current_node++); $$->set_sequence(scan.seq.get_taxa_index($1));} 
;

%%


void parse_tree(const char *file, const InputSequences &sequences, Tree &o_tree) {
  FILE *myfile = fopen(file, "r");
  if (!myfile) {
    cout << "Can't open " << file << endl;
    return;
  }
  rb_treein = myfile;
  TreeScanner scanner(sequences, o_tree);
  o_tree.init(sequences.number());
  do {
    rb_treeparse(scanner);
  } while (!feof(rb_treein));

}


void rb_treeerror(TreeScanner &scanner, const char *s) {
  cout << "Trees, parse error!  Message: " << &scanner << s << endl;
  // might as well halt now:
  exit(-1);
}
