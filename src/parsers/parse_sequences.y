%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "../repeatsbalance.h"
using namespace std;

// stuff from flex that bison needs to know about:
extern  int yylex();
extern  FILE *yyin;

void yyerror(InputSequences *oseq, const char *s);
%}


%union {
  int ival;
  float fval;
  char *sval;
  InputSequences *oseq;
}

%parse-param {InputSequences *oseq}

%token <ival> INT
%token <sval> STRING
%token NEWLINE
%type<oseq> sequences

%%



sequences:
   first_line  {;}
  | sequences STRING STRING NEWLINE { oseq->add_sequence($2, $3);}
  ;

first_line:
          INT INT NEWLINE {}
;

%%

void parse_sequences(const char *file, InputSequences &sequences) {
  // open a file handle to a particular file:
  FILE *myfile = fopen(file, "r");
  // make sure it is valid:
  if (!myfile) {
    cout << "Can't open " << file << endl;
    return;
  }
  // set flex to read from it instead of defaulting to STDIN:
  yyin = myfile;
  
  // parse through the input until there is no more:
  do {
    yyparse(&sequences);
  } while (!feof(yyin));

}


void yyerror(InputSequences *oseq, const char *s) {
  cout << "EEK, parse error!  Message: " << oseq << " " << s << endl;
  // might as well halt now:
  exit(-1);
}
