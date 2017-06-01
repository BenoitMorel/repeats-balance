%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "../repeatsbalance.h"  // to get the token types that we return
#include "parse_partitions.h"  // to get the token types that we return

%}

%option prefix="rb_part"
%option noyywrap
%option nounput
%option noinput
%%
[\-]           {return DASH;}
[,]           {return COMMA;}
[\n]           {return NEWLINE;}
[=]           {return EQUAL;}
[0-9]+            { rb_partlval.ival = atoi(rb_parttext); return INT; }
[a-zA-Z0-9_\.]+      {
  rb_partlval.sval = new char [ strlen(rb_parttext) + 1 ];
  strcpy(rb_partlval.sval, rb_parttext);
  return STRING;
}
.                ;
%%
