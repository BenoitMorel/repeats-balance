%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "../repeatsbalance.h"  // to get the token types that we return
#include "../Node.hpp"  // to get the token types that we return
#include "parse_tree.h"  // to get the token types that we return

%}

%option prefix="rb_tree"
%option noyywrap
%option nounput
%option noinput
%%
[(]           {return LEFT_PAR;}
[)]           {return RIGHT_PAR;}
[,]           {return COMMA;}
[;]           {return SEMI_COLON;}
[a-zA-Z0-9_]+      {
  rb_treelval.sval = new char [ strlen(rb_treetext) + 1 ];
  strcpy(rb_treelval.sval, rb_treetext);
  return STRING;
}
.                ;
%%
