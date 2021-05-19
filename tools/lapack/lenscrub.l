/* {definitions} */
iofun	"("[^;\{]*[;\{]
decl	"("[^)]*")"[,;]
any	[.]*
S	[ \t\n]*
cS	","{S}
len	[a-z][a-z0-9]*_len

%%
"s_stop"{decl}          |
"do_fio"{decl}          |
"s_cat"{iofun}          |
"s_copy"{iofun}         |
"s_stop"{iofun}         |
"s_cmp"{iofun}          |
"i_len"{iofun}          |
"len_trim__"{iofun}     |
"do_fio"{iofun}         |
"do_lio"{iofun}         { printf("%s", yytext); /* unchanged */ }
{any}"ilaenv_("         |
"dvout_("               |
"dmout_("               |
"ivout_("               |
"xerbla_("              |
[a-z]"tim"[a-z0-9]*"_(" |
[a-z]"prtb"[a-z0-9]"_(" {
                          register int c, paran_count = 1;
                          printf("%s", yytext); /* unchanged */
			  /* Loop until the correct closing paranthesis */
                          while (paran_count != 0) {
                              c = input();
                              if (c == '(') ++paran_count;
                              else if (c == ')') --paran_count;
                              putchar(c);
                          }
                        }
{cS}"("{S}ftnlen{S}")"{S}[1-9][0-9]* { ; /* omit -- f2c -A */ }
{cS}[1-9]([0-9])*L	{ ; /* omit */ }
{cS}ftnlen({S}{len})?	{ ; /* omit -- f2c -A */ }
^ftnlen" "{len}";\n"	{ ; /* omit -- f2c without -A or -C++ */ }
{cS}{len}		{ ; }
.			{ printf("%s", yytext); /* unchanged */ }
