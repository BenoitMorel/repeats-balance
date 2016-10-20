%{
#include <cstdio>
#include <iostream>
#include <vector>
#include "repeatsbalance.h"  // to get the token types that we return
#include "parse_sequences.h"  // to get the token types that we return

%}

%option noyywrap
%option nounput
%option noinput
%%
[\n]           {return NEWLINE;}
[0-9]+            { yylval.ival = atoi(yytext); return INT; }
[a-zA-Z0-9]+      {

  yylval.sval = new char [ strlen(yytext) + 1 ];
  strcpy(yylval.sval, yytext);
  return STRING;
}
.                ;
%%
