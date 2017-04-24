%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "../repeatsbalance.h"
using namespace std;

// stuff from flex that bison needs to know about:
extern  int rb_partlex();
extern  FILE *rb_partin;

void rb_parterror(InputPartitions *opart, const char *s);
%}


%union {
  int ival;
  float fval;
  char *sval;
  InputPartitions *opart;
}

%parse-param {InputPartitions *opart}

%token <ival> INT
%token <sval> STRING
%token COMMA
%token EQUAL
%token DASH
%token NEWLINE
%type<opart> partitions
%define api.prefix {rb_part}

%%



partitions:
   STRING COMMA STRING  EQUAL INT DASH INT NEWLINE { opart->add_partition($3, $5 - 1, $7 - $5 + 1); delete[] $3; delete[] $1;}
  | partitions STRING COMMA STRING  EQUAL INT DASH INT NEWLINE { opart->add_partition($4, $6 - 1, $8 - $6 + 1); delete [] $4; delete[] $2;}
  ;


%%

void parse_partitions(const char *file, InputPartitions &partitions) {
  // open a file handle to a particular file:
  FILE *myfile = fopen(file, "r");
  // make sure it is valid:
  if (!myfile) {
    cout << "Can't open " << file << endl;
    return;
  }
  // set flex to read from it instead of defaulting to STDIN:
  rb_partin = myfile;
  
  // parse through the input until there is no more:
  do {
    rb_partparse(&partitions);
  } while (!feof(rb_partin));
  fclose(myfile);
}


void rb_parterror(InputPartitions *opart, const char *s) {
  cout << "Partitions, parse error!  Message: " << opart << " " << s << endl;
  // might as well halt now:
  exit(-1);
}
