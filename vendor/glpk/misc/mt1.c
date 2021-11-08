/* mt1.c (0-1 knapsack problem; Martello & Toth algorithm) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  THIS CODE IS THE RESULT OF TRANSLATION OF THE FORTRAN SUBROUTINES
*  MT1 FROM THE BOOK:
*
*  SILVANO MARTELLO, PAOLO TOTH. KNAPSACK PROBLEMS: ALGORITHMS AND
*  COMPUTER IMPLEMENTATIONS. JOHN WILEY & SONS, 1990.
*
*  THE TRANSLATION HAS BEEN DONE WITH THE PERMISSION OF THE AUTHORS OF
*  THE ORIGINAL FORTRAN SUBROUTINES: SILVANO MARTELLO AND PAOLO TOTH.
*
*  The translation was made by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#line 1 ""
/*  -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#if 0 /* by mao */
#include "f2c.h"
#else
#include "env.h"
#include "mt1.h"

typedef int integer;
typedef float real;
#endif

#line 1 ""
/*<       SUBROUTINE MT1(N,P,W,C,Z,X,JDIM,JCK,XX,MIN,PSIGN,WSIGN,ZSIGN) >*/
#if 1 /* by mao */
static int chmt1_(int *, int *, int *, int *, int *, int *);

static
#endif
/* Subroutine */ int mt1_(integer *n, integer *p, integer *w, integer *c__,
	integer *z__, integer *x, integer *jdim, integer *jck, integer *xx,
	integer *min__, integer *psign, integer *wsign, integer *zsign)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real a, b;
    static integer j, r__, t, j1, n1, ch, ii, jj, kk, in, ll, ip, nn, iu, ii1,
	     chs, lim, lim1, diff, lold, mink;
    extern /* Subroutine */ int chmt1_(integer *, integer *, integer *,
	    integer *, integer *, integer *);
    static integer profit;


/* THIS SUBROUTINE SOLVES THE 0-1 SINGLE KNAPSACK PROBLEM */

/* MAXIMIZE  Z = P(1)*X(1) + ... + P(N)*X(N) */

/* SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) .LE. C , */
/*               X(J) = 0 OR 1  FOR J=1,...,N. */

/* THE PROGRAM IS INCLUDED IN THE VOLUME */
/*   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS */
/*   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990 */
/* AND IMPLEMENTS THE BRANCH-AND-BOUND ALGORITHM DESCRIBED IN */
/* SECTION  2.5.2 . */
/* THE PROGRAM DERIVES FROM AN EARLIER CODE PRESENTED IN */
/*  S. MARTELLO, P. TOTH, "ALGORITHM FOR THE SOLUTION OF THE 0-1 SINGLE */
/*  KNAPSACK PROBLEM", COMPUTING, 1978. */

/* THE INPUT PROBLEM MUST SATISFY THE CONDITIONS */

/*   1) 2 .LE. N .LE. JDIM - 1 ; */
/*   2) P(J), W(J), C  POSITIVE INTEGERS; */
/*   3) MAX (W(J)) .LE. C ; */
/*   4) W(1) + ... + W(N) .GT. C ; */
/*   5) P(J)/W(J) .GE. P(J+1)/W(J+1) FOR J=1,...,N-1. */

/* MT1 CALLS  1  PROCEDURE: CHMT1. */

/* THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS */
/* ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MT1. */
/* NO MACHINE-DEPENDENT CONSTANT IS USED. */
/* THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN */
/* AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE */
/* SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY). */
/* THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P. */
/* 9000/840. */

/* MT1 NEEDS  8  ARRAYS ( P ,  W ,  X ,  XX ,  MIN ,  PSIGN ,  WSIGN */
/*                        AND  ZSIGN ) OF LENGTH AT LEAST  N + 1 . */

/* MEANING OF THE INPUT PARAMETERS: */
/* N    = NUMBER OF ITEMS; */
/* P(J) = PROFIT OF ITEM  J  (J=1,...,N); */
/* W(J) = WEIGHT OF ITEM  J  (J=1,...,N); */
/* C    = CAPACITY OF THE KNAPSACK; */
/* JDIM = DIMENSION OF THE 8 ARRAYS; */
/* JCK  = 1 IF CHECK ON THE INPUT DATA IS DESIRED, */
/*      = 0 OTHERWISE. */

/* MEANING OF THE OUTPUT PARAMETERS: */
/* Z    = VALUE OF THE OPTIMAL SOLUTION IF  Z .GT. 0 , */
/*      = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI- */
/*        TION  - Z  IS VIOLATED; */
/* X(J) = 1 IF ITEM  J  IS IN THE OPTIMAL SOLUTION, */
/*      = 0 OTHERWISE. */

/* ARRAYS XX, MIN, PSIGN, WSIGN AND ZSIGN ARE DUMMY. */

/* ALL THE PARAMETERS ARE INTEGER. ON RETURN OF MT1 ALL THE INPUT */
/* PARAMETERS ARE UNCHANGED. */

/*<       INTEGER P(JDIM),W(JDIM),X(JDIM),C,Z >*/
/*<       INTEGER XX(JDIM),MIN(JDIM),PSIGN(JDIM),WSIGN(JDIM),ZSIGN(JDIM) >*/
/*<       INTEGER CH,CHS,DIFF,PROFIT,R,T >*/
/*<       Z = 0 >*/
#line 65 ""
    /* Parameter adjustments */
#line 65 ""
    --zsign;
#line 65 ""
    --wsign;
#line 65 ""
    --psign;
#line 65 ""
    --min__;
#line 65 ""
    --xx;
#line 65 ""
    --x;
#line 65 ""
    --w;
#line 65 ""
    --p;
#line 65 ""

#line 65 ""
    /* Function Body */
#line 65 ""
    *z__ = 0;
/*<       IF ( JCK .EQ. 1 ) CALL CHMT1(N,P,W,C,Z,JDIM) >*/
#line 66 ""
    if (*jck == 1) {
#line 66 ""
	chmt1_(n, &p[1], &w[1], c__, z__, jdim);
#line 66 ""
    }
/*<       IF ( Z .LT. 0 ) RETURN >*/
#line 67 ""
    if (*z__ < 0) {
#line 67 ""
	return 0;
#line 67 ""
    }
/* INITIALIZE. */
/*<       CH = C >*/
#line 69 ""
    ch = *c__;
/*<       IP = 0 >*/
#line 70 ""
    ip = 0;
/*<       CHS = CH >*/
#line 71 ""
    chs = ch;
/*<       DO 10 LL=1,N >*/
#line 72 ""
    i__1 = *n;
#line 72 ""
    for (ll = 1; ll <= i__1; ++ll) {
/*<         IF ( W(LL) .GT. CHS ) GO TO 20 >*/
#line 73 ""
	if (w[ll] > chs) {
#line 73 ""
	    goto L20;
#line 73 ""
	}
/*<         IP = IP + P(LL) >*/
#line 74 ""
	ip += p[ll];
/*<         CHS = CHS - W(LL) >*/
#line 75 ""
	chs -= w[ll];
/*<    10 CONTINUE >*/
#line 76 ""
/* L10: */
#line 76 ""
    }
/*<    20 LL = LL - 1 >*/
#line 77 ""
L20:
#line 77 ""
    --ll;
/*<       IF ( CHS .EQ. 0 ) GO TO 50 >*/
#line 78 ""
    if (chs == 0) {
#line 78 ""
	goto L50;
#line 78 ""
    }
/*<       P(N+1) = 0 >*/
#line 79 ""
    p[*n + 1] = 0;
/*<       W(N+1) = CH + 1 >*/
#line 80 ""
    w[*n + 1] = ch + 1;
/*<       LIM = IP + CHS*P(LL+2)/W(LL+2) >*/
#line 81 ""
    lim = ip + chs * p[ll + 2] / w[ll + 2];
/*<       A = W(LL+1) - CHS >*/
#line 82 ""
    a = (real) (w[ll + 1] - chs);
/*<       B = IP + P(LL+1) >*/
#line 83 ""
    b = (real) (ip + p[ll + 1]);
/*<       LIM1 = B - A*FLOAT(P(LL))/FLOAT(W(LL)) >*/
#line 84 ""
    lim1 = b - a * (real) p[ll] / (real) w[ll];
/*<       IF ( LIM1 .GT. LIM ) LIM = LIM1 >*/
#line 85 ""
    if (lim1 > lim) {
#line 85 ""
	lim = lim1;
#line 85 ""
    }
/*<       MINK = CH + 1 >*/
#line 86 ""
    mink = ch + 1;
/*<       MIN(N) = MINK >*/
#line 87 ""
    min__[*n] = mink;
/*<       DO 30 J=2,N >*/
#line 88 ""
    i__1 = *n;
#line 88 ""
    for (j = 2; j <= i__1; ++j) {
/*<         KK = N + 2 - J >*/
#line 89 ""
	kk = *n + 2 - j;
/*<         IF ( W(KK) .LT. MINK ) MINK = W(KK) >*/
#line 90 ""
	if (w[kk] < mink) {
#line 90 ""
	    mink = w[kk];
#line 90 ""
	}
/*<         MIN(KK-1) = MINK >*/
#line 91 ""
	min__[kk - 1] = mink;
/*<    30 CONTINUE >*/
#line 92 ""
/* L30: */
#line 92 ""
    }
/*<       DO 40 J=1,N >*/
#line 93 ""
    i__1 = *n;
#line 93 ""
    for (j = 1; j <= i__1; ++j) {
/*<         XX(J) = 0 >*/
#line 94 ""
	xx[j] = 0;
/*<    40 CONTINUE >*/
#line 95 ""
/* L40: */
#line 95 ""
    }
/*<       Z = 0 >*/
#line 96 ""
    *z__ = 0;
/*<       PROFIT = 0 >*/
#line 97 ""
    profit = 0;
/*<       LOLD = N >*/
#line 98 ""
    lold = *n;
/*<       II = 1 >*/
#line 99 ""
    ii = 1;
/*<       GO TO 170 >*/
#line 100 ""
    goto L170;
/*<    50 Z = IP >*/
#line 101 ""
L50:
#line 101 ""
    *z__ = ip;
/*<       DO 60 J=1,LL >*/
#line 102 ""
    i__1 = ll;
#line 102 ""
    for (j = 1; j <= i__1; ++j) {
/*<         X(J) = 1 >*/
#line 103 ""
	x[j] = 1;
/*<    60 CONTINUE >*/
#line 104 ""
/* L60: */
#line 104 ""
    }
/*<       NN = LL + 1 >*/
#line 105 ""
    nn = ll + 1;
/*<       DO 70 J=NN,N >*/
#line 106 ""
    i__1 = *n;
#line 106 ""
    for (j = nn; j <= i__1; ++j) {
/*<         X(J) = 0 >*/
#line 107 ""
	x[j] = 0;
/*<    70 CONTINUE >*/
#line 108 ""
/* L70: */
#line 108 ""
    }
/*<       RETURN >*/
#line 109 ""
    return 0;
/* TRY TO INSERT THE II-TH ITEM INTO THE CURRENT SOLUTION. */
/*<    80 IF ( W(II) .LE. CH ) GO TO 90 >*/
#line 111 ""
L80:
#line 111 ""
    if (w[ii] <= ch) {
#line 111 ""
	goto L90;
#line 111 ""
    }
/*<       II1 = II + 1 >*/
#line 112 ""
    ii1 = ii + 1;
/*<       IF ( Z .GE. CH*P(II1)/W(II1) + PROFIT ) GO TO 280 >*/
#line 113 ""
    if (*z__ >= ch * p[ii1] / w[ii1] + profit) {
#line 113 ""
	goto L280;
#line 113 ""
    }
/*<       II = II1 >*/
#line 114 ""
    ii = ii1;
/*<       GO TO 80 >*/
#line 115 ""
    goto L80;
/* BUILD A NEW CURRENT SOLUTION. */
/*<    90 IP = PSIGN(II) >*/
#line 117 ""
L90:
#line 117 ""
    ip = psign[ii];
/*<       CHS = CH - WSIGN(II) >*/
#line 118 ""
    chs = ch - wsign[ii];
/*<       IN = ZSIGN(II) >*/
#line 119 ""
    in = zsign[ii];
/*<       DO 100 LL=IN,N >*/
#line 120 ""
    i__1 = *n;
#line 120 ""
    for (ll = in; ll <= i__1; ++ll) {
/*<         IF ( W(LL) .GT. CHS ) GO TO 160 >*/
#line 121 ""
	if (w[ll] > chs) {
#line 121 ""
	    goto L160;
#line 121 ""
	}
/*<         IP = IP + P(LL) >*/
#line 122 ""
	ip += p[ll];
/*<         CHS = CHS - W(LL) >*/
#line 123 ""
	chs -= w[ll];
/*<   100 CONTINUE >*/
#line 124 ""
/* L100: */
#line 124 ""
    }
/*<       LL = N >*/
#line 125 ""
    ll = *n;
/*<   110 IF ( Z .GE. IP + PROFIT ) GO TO 280 >*/
#line 126 ""
L110:
#line 126 ""
    if (*z__ >= ip + profit) {
#line 126 ""
	goto L280;
#line 126 ""
    }
/*<       Z = IP + PROFIT >*/
#line 127 ""
    *z__ = ip + profit;
/*<       NN = II - 1 >*/
#line 128 ""
    nn = ii - 1;
/*<       DO 120 J=1,NN >*/
#line 129 ""
    i__1 = nn;
#line 129 ""
    for (j = 1; j <= i__1; ++j) {
/*<         X(J) = XX(J) >*/
#line 130 ""
	x[j] = xx[j];
/*<   120 CONTINUE >*/
#line 131 ""
/* L120: */
#line 131 ""
    }
/*<       DO 130 J=II,LL >*/
#line 132 ""
    i__1 = ll;
#line 132 ""
    for (j = ii; j <= i__1; ++j) {
/*<         X(J) = 1 >*/
#line 133 ""
	x[j] = 1;
/*<   130 CONTINUE >*/
#line 134 ""
/* L130: */
#line 134 ""
    }
/*<       IF ( LL .EQ. N ) GO TO 150 >*/
#line 135 ""
    if (ll == *n) {
#line 135 ""
	goto L150;
#line 135 ""
    }
/*<       NN = LL + 1 >*/
#line 136 ""
    nn = ll + 1;
/*<       DO 140 J=NN,N >*/
#line 137 ""
    i__1 = *n;
#line 137 ""
    for (j = nn; j <= i__1; ++j) {
/*<         X(J) = 0 >*/
#line 138 ""
	x[j] = 0;
/*<   140 CONTINUE >*/
#line 139 ""
/* L140: */
#line 139 ""
    }
/*<   150 IF ( Z .NE. LIM ) GO TO 280 >*/
#line 140 ""
L150:
#line 140 ""
    if (*z__ != lim) {
#line 140 ""
	goto L280;
#line 140 ""
    }
/*<       RETURN >*/
#line 141 ""
    return 0;
/*<   160 IU = CHS*P(LL)/W(LL) >*/
#line 142 ""
L160:
#line 142 ""
    iu = chs * p[ll] / w[ll];
/*<       LL = LL - 1 >*/
#line 143 ""
    --ll;
/*<       IF ( IU .EQ. 0 ) GO TO 110 >*/
#line 144 ""
    if (iu == 0) {
#line 144 ""
	goto L110;
#line 144 ""
    }
/*<       IF ( Z .GE. PROFIT + IP + IU ) GO TO 280 >*/
#line 145 ""
    if (*z__ >= profit + ip + iu) {
#line 145 ""
	goto L280;
#line 145 ""
    }
/* SAVE THE CURRENT SOLUTION. */
/*<   170 WSIGN(II) = CH - CHS >*/
#line 147 ""
L170:
#line 147 ""
    wsign[ii] = ch - chs;
/*<       PSIGN(II) = IP >*/
#line 148 ""
    psign[ii] = ip;
/*<       ZSIGN(II) = LL + 1 >*/
#line 149 ""
    zsign[ii] = ll + 1;
/*<       XX(II) = 1 >*/
#line 150 ""
    xx[ii] = 1;
/*<       NN = LL - 1 >*/
#line 151 ""
    nn = ll - 1;
/*<       IF ( NN .LT. II) GO TO 190 >*/
#line 152 ""
    if (nn < ii) {
#line 152 ""
	goto L190;
#line 152 ""
    }
/*<       DO 180 J=II,NN >*/
#line 153 ""
    i__1 = nn;
#line 153 ""
    for (j = ii; j <= i__1; ++j) {
/*<         WSIGN(J+1) = WSIGN(J) - W(J) >*/
#line 154 ""
	wsign[j + 1] = wsign[j] - w[j];
/*<         PSIGN(J+1) = PSIGN(J) - P(J) >*/
#line 155 ""
	psign[j + 1] = psign[j] - p[j];
/*<         ZSIGN(J+1) = LL + 1 >*/
#line 156 ""
	zsign[j + 1] = ll + 1;
/*<         XX(J+1) = 1 >*/
#line 157 ""
	xx[j + 1] = 1;
/*<   180 CONTINUE >*/
#line 158 ""
/* L180: */
#line 158 ""
    }
/*<   190 J1 = LL + 1 >*/
#line 159 ""
L190:
#line 159 ""
    j1 = ll + 1;
/*<       DO 200 J=J1,LOLD >*/
#line 160 ""
    i__1 = lold;
#line 160 ""
    for (j = j1; j <= i__1; ++j) {
/*<         WSIGN(J) = 0 >*/
#line 161 ""
	wsign[j] = 0;
/*<         PSIGN(J) = 0 >*/
#line 162 ""
	psign[j] = 0;
/*<         ZSIGN(J) = J >*/
#line 163 ""
	zsign[j] = j;
/*<   200 CONTINUE >*/
#line 164 ""
/* L200: */
#line 164 ""
    }
/*<       LOLD = LL >*/
#line 165 ""
    lold = ll;
/*<       CH = CHS >*/
#line 166 ""
    ch = chs;
/*<       PROFIT = PROFIT + IP >*/
#line 167 ""
    profit += ip;
/*<       IF ( LL - (N - 2) ) 240, 220, 210 >*/
#line 168 ""
    if ((i__1 = ll - (*n - 2)) < 0) {
#line 168 ""
	goto L240;
#line 168 ""
    } else if (i__1 == 0) {
#line 168 ""
	goto L220;
#line 168 ""
    } else {
#line 168 ""
	goto L210;
#line 168 ""
    }
/*<   210 II = N >*/
#line 169 ""
L210:
#line 169 ""
    ii = *n;
/*<       GO TO 250 >*/
#line 170 ""
    goto L250;
/*<   220 IF ( CH .LT. W(N) ) GO TO 230 >*/
#line 171 ""
L220:
#line 171 ""
    if (ch < w[*n]) {
#line 171 ""
	goto L230;
#line 171 ""
    }
/*<       CH = CH - W(N) >*/
#line 172 ""
    ch -= w[*n];
/*<       PROFIT = PROFIT + P(N) >*/
#line 173 ""
    profit += p[*n];
/*<       XX(N) = 1 >*/
#line 174 ""
    xx[*n] = 1;
/*<   230 II = N - 1 >*/
#line 175 ""
L230:
#line 175 ""
    ii = *n - 1;
/*<       GO TO 250 >*/
#line 176 ""
    goto L250;
/*<   240 II = LL + 2 >*/
#line 177 ""
L240:
#line 177 ""
    ii = ll + 2;
/*<       IF ( CH .GE. MIN(II-1) ) GO TO 80 >*/
#line 178 ""
    if (ch >= min__[ii - 1]) {
#line 178 ""
	goto L80;
#line 178 ""
    }
/* SAVE THE CURRENT OPTIMAL SOLUTION. */
/*<   250 IF ( Z .GE. PROFIT ) GO TO 270 >*/
#line 180 ""
L250:
#line 180 ""
    if (*z__ >= profit) {
#line 180 ""
	goto L270;
#line 180 ""
    }
/*<       Z = PROFIT >*/
#line 181 ""
    *z__ = profit;
/*<       DO 260 J=1,N >*/
#line 182 ""
    i__1 = *n;
#line 182 ""
    for (j = 1; j <= i__1; ++j) {
/*<         X(J) = XX(J) >*/
#line 183 ""
	x[j] = xx[j];
/*<   260 CONTINUE >*/
#line 184 ""
/* L260: */
#line 184 ""
    }
/*<       IF ( Z .EQ. LIM ) RETURN >*/
#line 185 ""
    if (*z__ == lim) {
#line 185 ""
	return 0;
#line 185 ""
    }
/*<   270 IF ( XX(N) .EQ. 0 ) GO TO 280 >*/
#line 186 ""
L270:
#line 186 ""
    if (xx[*n] == 0) {
#line 186 ""
	goto L280;
#line 186 ""
    }
/*<       XX(N) = 0 >*/
#line 187 ""
    xx[*n] = 0;
/*<       CH = CH + W(N) >*/
#line 188 ""
    ch += w[*n];
/*<       PROFIT = PROFIT - P(N) >*/
#line 189 ""
    profit -= p[*n];
/* BACKTRACK. */
/*<   280 NN = II - 1 >*/
#line 191 ""
L280:
#line 191 ""
    nn = ii - 1;
/*<       IF ( NN .EQ. 0 ) RETURN >*/
#line 192 ""
    if (nn == 0) {
#line 192 ""
	return 0;
#line 192 ""
    }
/*<       DO 290 J=1,NN >*/
#line 193 ""
    i__1 = nn;
#line 193 ""
    for (j = 1; j <= i__1; ++j) {
/*<         KK = II - J >*/
#line 194 ""
	kk = ii - j;
/*<         IF ( XX(KK) .EQ. 1 ) GO TO 300 >*/
#line 195 ""
	if (xx[kk] == 1) {
#line 195 ""
	    goto L300;
#line 195 ""
	}
/*<   290 CONTINUE >*/
#line 196 ""
/* L290: */
#line 196 ""
    }
/*<       RETURN >*/
#line 197 ""
    return 0;
/*<   300 R = CH >*/
#line 198 ""
L300:
#line 198 ""
    r__ = ch;
/*<       CH = CH + W(KK) >*/
#line 199 ""
    ch += w[kk];
/*<       PROFIT = PROFIT - P(KK) >*/
#line 200 ""
    profit -= p[kk];
/*<       XX(KK) = 0 >*/
#line 201 ""
    xx[kk] = 0;
/*<       IF ( R .LT. MIN(KK) ) GO TO 310 >*/
#line 202 ""
    if (r__ < min__[kk]) {
#line 202 ""
	goto L310;
#line 202 ""
    }
/*<       II = KK + 1 >*/
#line 203 ""
    ii = kk + 1;
/*<       GO TO 80 >*/
#line 204 ""
    goto L80;
/*<   310 NN = KK + 1 >*/
#line 205 ""
L310:
#line 205 ""
    nn = kk + 1;
/*<       II = KK >*/
#line 206 ""
    ii = kk;
/* TRY TO SUBSTITUTE THE NN-TH ITEM FOR THE KK-TH. */
/*<   320 IF ( Z .GE. PROFIT + CH*P(NN)/W(NN) ) GO TO 280 >*/
#line 208 ""
L320:
#line 208 ""
    if (*z__ >= profit + ch * p[nn] / w[nn]) {
#line 208 ""
	goto L280;
#line 208 ""
    }
/*<       DIFF = W(NN) - W(KK) >*/
#line 209 ""
    diff = w[nn] - w[kk];
/*<       IF ( DIFF ) 370, 330, 340 >*/
#line 210 ""
    if (diff < 0) {
#line 210 ""
	goto L370;
#line 210 ""
    } else if (diff == 0) {
#line 210 ""
	goto L330;
#line 210 ""
    } else {
#line 210 ""
	goto L340;
#line 210 ""
    }
/*<   330 NN = NN + 1 >*/
#line 211 ""
L330:
#line 211 ""
    ++nn;
/*<       GO TO 320 >*/
#line 212 ""
    goto L320;
/*<   340 IF ( DIFF .GT. R ) GO TO 330 >*/
#line 213 ""
L340:
#line 213 ""
    if (diff > r__) {
#line 213 ""
	goto L330;
#line 213 ""
    }
/*<       IF ( Z .GE. PROFIT + P(NN) ) GO TO 330 >*/
#line 214 ""
    if (*z__ >= profit + p[nn]) {
#line 214 ""
	goto L330;
#line 214 ""
    }
/*<       Z = PROFIT + P(NN) >*/
#line 215 ""
    *z__ = profit + p[nn];
/*<       DO 350 J=1,KK >*/
#line 216 ""
    i__1 = kk;
#line 216 ""
    for (j = 1; j <= i__1; ++j) {
/*<         X(J) = XX(J) >*/
#line 217 ""
	x[j] = xx[j];
/*<   350 CONTINUE >*/
#line 218 ""
/* L350: */
#line 218 ""
    }
/*<       JJ = KK + 1 >*/
#line 219 ""
    jj = kk + 1;
/*<       DO 360 J=JJ,N >*/
#line 220 ""
    i__1 = *n;
#line 220 ""
    for (j = jj; j <= i__1; ++j) {
/*<         X(J) = 0 >*/
#line 221 ""
	x[j] = 0;
/*<   360 CONTINUE >*/
#line 222 ""
/* L360: */
#line 222 ""
    }
/*<       X(NN) = 1 >*/
#line 223 ""
    x[nn] = 1;
/*<       IF ( Z .EQ. LIM ) RETURN >*/
#line 224 ""
    if (*z__ == lim) {
#line 224 ""
	return 0;
#line 224 ""
    }
/*<       R = R - DIFF >*/
#line 225 ""
    r__ -= diff;
/*<       KK = NN >*/
#line 226 ""
    kk = nn;
/*<       NN = NN + 1 >*/
#line 227 ""
    ++nn;
/*<       GO TO 320 >*/
#line 228 ""
    goto L320;
/*<   370 T = R - DIFF >*/
#line 229 ""
L370:
#line 229 ""
    t = r__ - diff;
/*<       IF ( T .LT. MIN(NN) ) GO TO 330 >*/
#line 230 ""
    if (t < min__[nn]) {
#line 230 ""
	goto L330;
#line 230 ""
    }
/*<       IF ( Z .GE. PROFIT + P(NN) + T*P(NN+1)/W(NN+1)) GO TO 280 >*/
#line 231 ""
    if (*z__ >= profit + p[nn] + t * p[nn + 1] / w[nn + 1]) {
#line 231 ""
	goto L280;
#line 231 ""
    }
/*<       CH = CH - W(NN) >*/
#line 232 ""
    ch -= w[nn];
/*<       PROFIT = PROFIT + P(NN) >*/
#line 233 ""
    profit += p[nn];
/*<       XX(NN) = 1 >*/
#line 234 ""
    xx[nn] = 1;
/*<       II = NN + 1 >*/
#line 235 ""
    ii = nn + 1;
/*<       WSIGN(NN) = W(NN) >*/
#line 236 ""
    wsign[nn] = w[nn];
/*<       PSIGN(NN) = P(NN) >*/
#line 237 ""
    psign[nn] = p[nn];
/*<       ZSIGN(NN) = II >*/
#line 238 ""
    zsign[nn] = ii;
/*<       N1 = NN + 1 >*/
#line 239 ""
    n1 = nn + 1;
/*<       DO 380 J=N1,LOLD >*/
#line 240 ""
    i__1 = lold;
#line 240 ""
    for (j = n1; j <= i__1; ++j) {
/*<         WSIGN(J) = 0 >*/
#line 241 ""
	wsign[j] = 0;
/*<         PSIGN(J) = 0 >*/
#line 242 ""
	psign[j] = 0;
/*<         ZSIGN(J) = J >*/
#line 243 ""
	zsign[j] = j;
/*<   380 CONTINUE >*/
#line 244 ""
/* L380: */
#line 244 ""
    }
/*<       LOLD = NN >*/
#line 245 ""
    lold = nn;
/*<       GO TO 80 >*/
#line 246 ""
    goto L80;
/*<       END >*/
} /* mt1_ */

/*<       SUBROUTINE CHMT1(N,P,W,C,Z,JDIM) >*/
#if 1 /* by mao */
static
#endif
/* Subroutine */ int chmt1_(integer *n, integer *p, integer *w, integer *c__,
	integer *z__, integer *jdim)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static real r__, rr;
    static integer jsw;


/* CHECK THE INPUT DATA. */

/*<       INTEGER P(JDIM),W(JDIM),C,Z >*/
/*<       IF ( N .GE. 2 .AND. N .LE. JDIM - 1 ) GO TO 10 >*/
#line 253 ""
    /* Parameter adjustments */
#line 253 ""
    --w;
#line 253 ""
    --p;
#line 253 ""

#line 253 ""
    /* Function Body */
#line 253 ""
    if (*n >= 2 && *n <= *jdim - 1) {
#line 253 ""
	goto L10;
#line 253 ""
    }
/*<       Z = - 1 >*/
#line 254 ""
    *z__ = -1;
/*<       RETURN >*/
#line 255 ""
    return 0;
/*<    10 IF ( C .GT. 0 ) GO TO 30 >*/
#line 256 ""
L10:
#line 256 ""
    if (*c__ > 0) {
#line 256 ""
	goto L30;
#line 256 ""
    }
/*<    20 Z = - 2 >*/
#line 257 ""
L20:
#line 257 ""
    *z__ = -2;
/*<       RETURN >*/
#line 258 ""
    return 0;
/*<    30 JSW = 0 >*/
#line 259 ""
L30:
#line 259 ""
    jsw = 0;
/*<       RR = FLOAT(P(1))/FLOAT(W(1)) >*/
#line 260 ""
    rr = (real) p[1] / (real) w[1];
/*<       DO 50 J=1,N >*/
#line 261 ""
    i__1 = *n;
#line 261 ""
    for (j = 1; j <= i__1; ++j) {
/*<         R = RR >*/
#line 262 ""
	r__ = rr;
/*<         IF ( P(J) .LE. 0 ) GO TO 20 >*/
#line 263 ""
	if (p[j] <= 0) {
#line 263 ""
	    goto L20;
#line 263 ""
	}
/*<         IF ( W(J) .LE. 0 ) GO TO 20 >*/
#line 264 ""
	if (w[j] <= 0) {
#line 264 ""
	    goto L20;
#line 264 ""
	}
/*<         JSW = JSW + W(J) >*/
#line 265 ""
	jsw += w[j];
/*<         IF ( W(J) .LE. C ) GO TO 40 >*/
#line 266 ""
	if (w[j] <= *c__) {
#line 266 ""
	    goto L40;
#line 266 ""
	}
/*<         Z = - 3 >*/
#line 267 ""
	*z__ = -3;
/*<         RETURN >*/
#line 268 ""
	return 0;
/*<    40   RR = FLOAT(P(J))/FLOAT(W(J)) >*/
#line 269 ""
L40:
#line 269 ""
	rr = (real) p[j] / (real) w[j];
/*<         IF ( RR .LE. R ) GO TO 50 >*/
#line 270 ""
	if (rr <= r__) {
#line 270 ""
	    goto L50;
#line 270 ""
	}
/*<         Z = - 5 >*/
#line 271 ""
	*z__ = -5;
/*<         RETURN >*/
#line 272 ""
	return 0;
/*<    50 CONTINUE >*/
#line 273 ""
L50:
#line 273 ""
	;
#line 273 ""
    }
/*<       IF ( JSW .GT. C ) RETURN >*/
#line 274 ""
    if (jsw > *c__) {
#line 274 ""
	return 0;
#line 274 ""
    }
/*<       Z = - 4 >*/
#line 275 ""
    *z__ = -4;
/*<       RETURN >*/
#line 276 ""
    return 0;
/*<       END >*/
} /* chmt1_ */

#if 1 /* by mao */
int mt1(int n, int p[], int w[], int c, int x[], int jck, int xx[],
      int min[], int psign[], int wsign[], int zsign[])
{     /* solve 0-1 knapsack problem */
      int z, jdim = n+1, j, s1, s2;
      mt1_(&n, &p[1], &w[1], &c, &z, &x[1], &jdim, &jck, &xx[1],
         &min[1], &psign[1], &wsign[1], &zsign[1]);
      /* check solution found */
      s1 = s2 = 0;
      for (j = 1; j <= n; j++)
      {  xassert(x[j] == 0 || x[j] == 1);
         if (x[j])
            s1 += p[j], s2 += w[j];
      }
      xassert(s1 == z);
      xassert(s2 <= c);
      return z;
}
#endif

/* eof */
