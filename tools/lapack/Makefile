LOADLIBS = -ly -lfl -lm  
LIBS = -lfl -lm  
CFLAGS	= -O

all: lenscrub comment

lenscrub: lenscrub.l
	lex lenscrub.l
	mv -f lex.yy.c lex_for_lenscrub.c
	cc -o lenscrub -O lex_for_lenscrub.c -ll

comment: comment.l
	lex comment.l
	mv -f lex.yy.c lex_for_comment.c
	cc -o comment -O lex_for_comment.c -ll
