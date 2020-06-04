/******************************************************************************
 * bn.c - big number math implementation
 *
 * Copyright (c) 2004 by Juergen Buchmueller <pullmoll@stop1984.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
 *  $Id: bignum.c,v 1.17 2005/07/23 02:55:53 pullmoll Exp $
 ******************************************************************************/

#include "bignum.h"
#include "igraph_error.h"
#include "config.h"

#ifndef ASM_X86
    #ifdef  X86
        #define ASM_X86 1
    #endif
#endif

/**
 * @brief Return hex representation of a big number
 *
 * Returns the hex representation of a[],
 * where a is a big number integer with nlimb limbs.
 *
 * @param a pointer to an array of limbs
 * @param nlimb number of limbs in the array
 *
 * @result string containing the hex representation of a
 */
const char *bn2x(limb_t *a, count_t nlimb) {
    static IGRAPH_THREAD_LOCAL count_t which = 0;
    static IGRAPH_THREAD_LOCAL char *xbuff[8] = {
        NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL
    };
    char *dst;
    count_t size;
    count_t n = nlimb;

    if (0 == n) {
        return "0";
    }

    which = (which + 1) % 8;
    size = 8 * n + 1;
    if (NULL != xbuff[which]) {
        free(xbuff[which]);
    }
    dst = xbuff[which] = calloc(size, sizeof(char));
    if (NULL == dst) {
        return "memory error";
    }
    while (n-- > 0) {
        dst += snprintf(dst, size, "%08x", a[n]);
        size -= 8;
    }
    return xbuff[which];
}

/**
 * @brief Return decimal representation of a big number
 *
 * Returns the decimal representation of a[],
 * where a is a big number integer with nlimb limbs.
 *
 * @param a pointer to an array of limbs
 * @param nlimb number of limbs in the array
 *
 * @result string containing the decimal representation of a
 */
const char *bn2d(limb_t *a, count_t nlimb) {
    static IGRAPH_THREAD_LOCAL count_t which = 0;
    static IGRAPH_THREAD_LOCAL char *dbuff[8] = {
        NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL
    };
    static IGRAPH_THREAD_LOCAL limb_t v[BN_MAXSIZE];
    limb_t r;
    char *dst;
    count_t size;
    count_t n = bn_sizeof(a, nlimb);

    if (0 == n) {
        return "0";
    }

    bn_copy(v, a, n);
    which = (which + 1) % 8;
    size = 12 * n + 1;
    if (NULL != dbuff[which]) {
        free(dbuff[which]);
    }
    dst = dbuff[which] = calloc(size, sizeof(char));
    if (NULL == dst) {
        return "memory error";
    }
    size--;
    while (0 != bn_cmp_limb(v, 0, n)) {
        r = bn_div_limb(v, v, 10, n);
        dst[--size] = '0' + (char) r;
    }
    return &dst[size];
}

/**
 * @brief Return decimal representation of a big number pair
 *
 * Returns the decimal representation of a[].b[],
 * where a is a big number integer with alimb limbs,
 * and b is a multiprecision fixed fraction with blimb limbs.
 *
 * @param a pointer to an array of limbs
 * @param alimb number of limbs in the a array
 * @param b pointer to an array of limbs
 * @param blimb number of limbs in the b array
 *
 * @result string containing the decimal representation of a.b
 */
const char *bn2f(limb_t *a, count_t alimb, limb_t *b, count_t blimb) {
    static IGRAPH_THREAD_LOCAL count_t which = 0;
    static IGRAPH_THREAD_LOCAL char *dbuff[8] = {
        NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL
    };
    static IGRAPH_THREAD_LOCAL limb_t v[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t w[BN_MAXSIZE];
    limb_t r;
    char *dst;
    count_t size;

    bn_copy(v, a, alimb);
    bn_copy(w, b, blimb);

    which = (which + 1) % 8;
    size = 12 * (alimb + blimb) + 1 + 1;
    if (NULL != dbuff[which]) {
        free(dbuff[which]);
    }
    dst = dbuff[which] = calloc(size, sizeof(char));
    if (NULL == dst) {
        return "memory error";
    }
    size = 12 * alimb;
    while (0 != bn_cmp_limb(w, 0, blimb) && size < 12 * (alimb + blimb)) {
        r = bn_mul_limb(w, w, 10, blimb);
        dst[size++] = '0' + (char) r;
    }

    size = 12 * alimb;
    dst[size] = '.';
    while (0 != bn_cmp_limb(v, 0, alimb) && size > 0) {
        r = bn_div_limb(v, v, 10, alimb);
        dst[--size] = '0' + (char) r;
    }

    return &dst[size];
}

/**
 * @brief Return binary representation of a big number
 *
 * Returns the binary representation of a[],
 * where a is a big number integer with nlimb limbs.
 *
 * @param a pointer to an array of limbs
 * @param nlimb number of limbs in the array
 *
 * @result string containing the binary representation of a
 */
const char *bn2b(limb_t *a, count_t nlimb) {
    static IGRAPH_THREAD_LOCAL count_t which = 0;
    static IGRAPH_THREAD_LOCAL char *bbuff[8] = {
        NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL
    };
    limb_t r;
    char *dst;
    count_t size;
    count_t n = bn_sizeof(a, nlimb);

    if (0 == n) {
        return "0";
    }

    which = (which + 1) % 8;
    size = LIMBBITS * n + 1;
    if (NULL != bbuff[which]) {
        free(bbuff[which]);
    }
    dst = bbuff[which] = calloc(size, sizeof(char));
    if (NULL == dst) {
        return "memory error";
    }
    n = 0;
    size--;
    while (size-- > 0) {
        r = (a[n / LIMBBITS] >> (n % LIMBBITS)) & 1;
        n++;
        dst[size] = '0' + (char) r;
    }
    return &dst[size];
}

/**
 * @brief Zero an array of limbs
 *
 * Sets a[] = 0
 * where a is a big number integer of nlimb limbs.
 *
 * @param a pointer to an array of limbs
 * @param nlimb number of limbs in the array
 *
 */
void bn_zero(limb_t a[], count_t nlimb) {
    memset(a, 0, nlimb * sizeof(limb_t));
}

/**
 * @brief Set an array of limbs to a single limb value
 *
 * Sets a[] = d
 * where a is a big number integer of nlimb limbs,
 * and d is a single limb
 *
 * @param a pointer to an array of limbs to set
 * @param d limb value to set a to
 * @param nlimb number of limbs in the array
 *
 */
void bn_limb(limb_t a[], limb_t d, count_t nlimb) {
    memset(a, 0, nlimb * sizeof(limb_t));
    a[0] = d;
}

/**
 * @brief Copy an array of limbs
 *
 * Sets a[] = b[]
 * where a and b are a big number integers of nlimb limbs
 *
 * @param a pointer to an array of limbs (destination)
 * @param b pointer to an array of limbs (source)
 * @param nlimb number of limbs in the arrays
 */
void bn_copy(limb_t a[], limb_t b[], count_t nlimb) {
    memcpy(a, b, nlimb * sizeof(limb_t));
}

/**
 * @brief Return significant size of a big number
 *
 * Returns size of significant limbs in a[]
 * i.e. searches for the first non-zero limb from
 * nlimb-1 downto 0.
 *
 * @param a pointer to an array of limbs (candidate)
 * @param nlimb number of limbs in the arrays
 *
 * @result number of significant limbs in a
 */
count_t bn_sizeof(limb_t a[], count_t nlimb) {
    while (nlimb-- > 0)
        if (0 != a[nlimb]) {
            return ++nlimb;
        }
    return 0;
}


/**
 * @brief Return sign of a bignum minus a limb
 *
 * Returns the sign of (a[] - b)
 * where a is a big number integer of nlimb limbs,
 * and b is a single limb
 +
 * @param a pointer to an array of limbs (minuend)
 * @param b a single limb (subtrahend)
 * @param nlimb number of limbs in the array a
 *
 * @result sign of the comparison: -1 a<b, 0 a=b, +1 a>b
 */
int bn_cmp_limb(limb_t a[], limb_t b, count_t nlimb) {
    if (0 == nlimb) {
        return 0;
    }

    while (nlimb-- > 1)
        if (0 != a[nlimb]) {
            return +1;
        }
    if (a[0] < b) {
        return -1;
    }
    if (a[0] > b) {
        return +1;
    }
    return 0;
}

/**
 * @brief Return sign of bignum a minus bignum b
 *
 * Returns the sign of (a[] - b[])
 * where a and b are a big number integers of nlimb limbs
 *
 * @param a pointer to an array of limbs (minuend)
 * @param b pointer to an array of limbs (subtrahend)
 * @param nlimb number of limbs in the arrays
 *
 * @result sign of the comparison: -1 a<b, 0 a=b, +1 a>b
 */
int bn_cmp(limb_t a[], limb_t b[], count_t nlimb) {
    if (0 == nlimb) {
        return 0;
    }

    while (nlimb-- > 0) {
        if (a[nlimb] > b[nlimb]) {
            return +1;    /* GT */
        }
        if (a[nlimb] < b[nlimb]) {
            return -1;    /* LT */
        }
    }

    return 0;   /* EQ */
}

/**
 * @brief Single limb is even test
 *
 * Returns 1 if a is even, else 0
 * where a is a single limb
 *
 * @param a a single limb
 *
 * @result zero if a is odd, 1 if a is even
 */
int sl_iseven(limb_t a) {
    return (a & 1) ? 0 : 1;
}

/**
 * @brief bignum is even test
 *
 * Returns 1 if a[] is even, else 0
 * where a is a big number integer of nlimb limbs
 * Note: a zero limb big number integer is even!
 *
 * @param a pointer to an array of limbs
 * @param nlimb number of limbs in the arrays
 *
 * @result zero if a is odd, 1 if a is even
 */
int bn_iseven(limb_t *a, count_t nlimb) {
    if (0 == nlimb) {
        return 1;
    }
    return (a[0] & 1) ? 0 : 1;
}

/**
 * @brief Add a single limb to a bignum
 *
 * Computes w[] = u[] + v
 * where w, u are big number integers of nlimb lims each,
 * and v is a single limb.
 * Returns carry if the addition overflows.
 *
 * Ref: Derived from Knuth Algorithm A.
 *
 * @param w pointer to an array of limbs receiving result
 * @param u pointer to an array of limbs (addend 1)
 * @param v a single limb
 * @param nlimb number of limbs in the arrays w and u
 *
 * @result The carry status of the addition
 */
limb_t bn_add_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb) {
    limb_t carry;
    count_t j;

    /* Copy u to w, so we can bail out if no borrow is left */
    if (w != u) {
        bn_copy(w, u, nlimb);
    }

    /* Add v to first limb of u */
    w[0] += v;
    carry = (w[0] < v ? 1 : 0);

    /* Add carry to subsequent limbs */
    for (j = 1; 0 != carry && j < nlimb; j++) {
        w[j] += carry;
        carry = (w[j] < carry ? 1 : 0);
    }
    return carry;
}


/**
 * @brief Subtract a single limb from a bignum
 *
 * Computes w[] = u[] - v
 * where w, u are big number integers of nlimb limbs each,
 * and v is a single limb.
 * Returns borrow (0 if u >= v, or 1 if v > u).
 *
 * Ref: Derived from Knuth Algorithm S.
 *
 * @param w pointer to an array of limbs receiving the result
 * @param u pointer to an array of limbs (minuend)
 * @param v single limb (subtrahend)
 * @param nlimb number of limbs in the arrays
 *
 * @result borrow of the subtraction (0 if u >= v, 1 if u < v)
 */
limb_t bn_sub_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb) {
    limb_t borrow;
    count_t j;

    /* Copy u to w, so we can bail out if no borrow is left */
    if (w != u) {
        bn_copy(w, u, nlimb);
    }

    /* Subtract v from first limb of u */
    w[0] -= v;
    borrow = (w[0] > ~v ? 1 : 0);

    /* Subtract borrow from subsequent limbs */
    for (j = 1; 0 != borrow && j < nlimb; j++) {
        w[j] -= borrow;
        borrow = (w[j] > ~borrow ? 1 : 0);
    }

    return borrow;
}

/**
 * @brief Divide a bignum by a single limb
 *
 * Computes quotient q[] = u[] / v
 * and returns remainder r = u[] % v
 * where q, u are big number integers of nlimb limbs each,
 * and v is a single limb.
 *
 * Makes no assumptions about normalisation.
 *
 * Ref: Knuth Vol 2 Ch 4.3.1 Exercise 16 p625
 *
 * @param q pointer to an array of limbs receiving the quotient
 * @param u pointer to an array of limbs (dividend)
 * @param v single limb (divisor)
 * @param nlimb number of limbs in the arrays
 *
 * @result single limb remainder of the division (modulo)
 */
limb_t bn_div_limb(limb_t q[], limb_t u[], limb_t v, count_t nlimb) {
    count_t j;
    limb_t t[2], r;
    count_t shift;

    if (0 == nlimb) {
        return 0;
    }
    if (0 == v) {
        return LIMBMASK;    /* Divide by zero error */
    }

    /*
     * Normalize first:
     * qequires high bit of V to be set,
     * so find most significant by shifting
     * until DIGMSB is set.
     */
    for (shift = 0; 0 == (v & DIGMSB); shift++) {
        v <<= 1;
    }
    r = bn_shl(q, u, shift, nlimb);

    j = nlimb;
    while (j-- > 0) {
        t[0] = q[j];
        t[1] = r;
        sl_div(&q[j], &r, t, v);
    }

    /* Unnormalize */
    r >>= shift;
    return r;
}

/**
 * @brief Modulo a bignum by a single limb
 *
 * Computes remainder (modulo) r = u[] mod v
 * Computes r = u[] mod v
 * where u is a big number integer of nlimb
 * and r, v are single precision limbs
 *
 * Use remainder from divide function.
 *
 * @param u pointer to an array of limbs (dividend)
 * @param v single limb (divisor)
 * @param nlimb number of limbs in the arrays
 *
 * @result single limb remainder of the division (modulo)
 */
limb_t bn_mod_limb(limb_t u[], limb_t v, count_t nlimb) {
    static IGRAPH_THREAD_LOCAL limb_t q[2 * BN_MAXSIZE];
    limb_t r;

    r = bn_div_limb(q, u, v, nlimb);

    bn_zero(q, nlimb);
    return r;
}

/**
 * @brief Multiply a bignum by a single limb
 *
 * Computes product w[] = u[] * v
 * Returns overflow k
 * where w, u are big number integers of nlimb each
 * and v is a single limb
 *
 * @param w pointer to an array of limbs to receive the result
 * @param u pointer to an array of limbs (factor)
 * @param v single limb (other factor)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero if no overflow, else overflow (value of w[nlimb])
 */
limb_t bn_mul_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb) {
    limb_t t[2];
    limb_t carry;
    count_t j;

    if (0 == v) {
        bn_zero(w, nlimb);
        return 0;
    }

    for (j = 0, carry = 0; j < nlimb; j++) {
        sl_mul(t, u[j], v);
        w[j] = t[0] + carry;
        carry = t[1] + (w[j] < carry ? 1 : 0);
    }

    return carry;
}

#if HAVE_U64
/**
 * @brief Computes quotient and remainder of 64 bit / 32 bit
 *
 * Computes quotient q = u[] / v, remainder r = u[] mod v
 * where u[] is a double limb.
 *
 * With native support for double limb division
 *
 * @param q pointer to the limb to receive the quotient
 * @param r pointer to the limb to receive the remainder
 * @param u pointer to an array of two limbs
 * @param v single limb divisor
 *
 * @result zero on success
 */
limb_t sl_div(limb_t *q, limb_t *r, limb_t u[2], limb_t v) {
#if ASM_X86
    limb_t qq;
    limb_t rr;

    if (0 == v)
        /* division by zero */
    {
        return LIMBMASK;
    }
    asm volatile(
        "divl	%4"
        : "=a"(qq), "=d"(rr)
        : "a"(u[0]), "d"(u[1]), "g"(v));
    *q = qq;
    *r = rr;
#else
    dlimb_t dd;

    if (0 == v)
        /* division by zero */
    {
        return LIMBMASK;
    }
    dd = ((dlimb_t)u[1] << LIMBBITS) | u[0];
    *q = (limb_t) (dd / v);
    *r = dd % v;
#endif
    return 0;
}

#else

#define B (HALFMASK + 1)

/**
 * @brief Computes quotient and remainder of 64 bit / 32 bit
 *
 * Computes quotient q = u / v, remainder r = u mod v
 * where u is a double limb
 * and q, v, r are single precision limbs.
 * Returns high limb of quotient (max value is 1)
 * Assumes normalized such that v1 >= b/2
 * where b is size of HALF_DIGIT
 * i.e. the most significant bit of v should be one
 *
 * In terms of half-limbs in Knuth notation:
 *   (q2q1q0) = (u4u3u2u1u0) / (v1v0)
 *   (r1r0) = (u4u3u2u1u0) % (v1v0)
 * for m = 2, n = 2 where u4 = 0
 *
 * We set q = (q1q0) and return q2 as "overflow'
 * Returned q2 is either 0 or 1.
 *
 * @param q pointer to the limb to receive the quotient
 * @param r pointer to the limb to receive the remainder
 * @param u pointer to an array of two limbs
 * @param v single limb divisor
 *
 * @result zero on success
 */
limb_t sl_div(limb_t *q, limb_t *r, limb_t u[2], limb_t v) {
    limb_t quot;
    limb_t rem;
    limb_t ul;
    limb_t uh;
    limb_t p0;
    limb_t p1;
    limb_t v0;
    limb_t v1;
    limb_t u0;
    limb_t u1;
    limb_t u2;
    limb_t u3;
    limb_t borrow;
    limb_t q1;
    limb_t q2;
    limb_t s;
    limb_t t;

    /* Check for normalisation */
    if (0 == (v & DIGMSB)) {
        *q = *r = 0;
        return LIMBMASK;
    }

    /* Split up into half-limbs */
    v0 = LSH(v);
    v1 = MSH(v);
    u0 = LSH(u[0]);
    u1 = MSH(u[0]);
    u2 = LSH(u[1]);
    u3 = MSH(u[1]);

    /* Do three rounds of Knuth Algorithm D Vol 2 p272 */

    /*
     * ROUND 1 calculate q2:
     * estimate quot = (u4u3)/v1 = 0 or 1,
     * then set (u4u3u2) -= quot*(v1v0) where u4 = 0.
     */
    quot = u3 / v1;
    if (quot > 0) {
        rem = u3 - quot * v1;
        t = SHL(rem) | u2;
        if (quot * v0 > t) {
            quot--;
        }
    }
    uh = 0;     /* (u4) */
    ul = u[1];  /* (u3u2) */
    if (quot > 0) {
        /* (u4u3u2) -= quot*(v1v0) where u4 = 0 */
        p0 = quot * v0;
        p1 = quot * v1;
        s = p0 + SHL(p1);
        ul -= s;
        borrow = (ul > ~s ? 1 : 0);
        uh -= MSH(p1) - borrow;

        if (0 != MSH(uh)) {
            /* add back */
            quot--;
            ul += v;
            uh = 0;
        }
    }
    q2 = quot;

    /*
     * ROUND 2 calculate q1:
     * estimate quot = (u3u2) / v1,
     * then set (u3u2u1) -= quot*(v1v0)
     */
    t = ul;
    quot = t / v1;
    rem = t - quot * v1;
    /* Test on v0 */
    t = SHL(rem) | u1;
    if (B == quot || (quot * v0) > t) {
        quot--;
        rem += v1;
        t = SHL(rem) | u1;
        if (rem < B && (quot * v0) > t) {
            quot--;
        }
    }

    /*
     * multiply and subtract:
     * (u3u2u1)' = (u3u2u1) - quot*(v1v0)
     */
    uh = MSH(ul);   /* (0u3) */
    ul = SHL(ul) | u1;  /* (u2u1) */
    p0 = quot * v0;
    p1 = quot * v1;
    s = p0 + SHL(p1);
    ul -= s;
    borrow = (ul > ~s ? 1 : 0);
    uh -= MSH(p1) - borrow;

    if (0 != MSH(uh)) {
        /* add back v */
        quot--;
        ul += v;
        uh = 0;
    }

    /* quotient q1 */
    q1 = quot;

    /*
     * ROUND 3:
     * calculate q0; estimate quot = (u2u1) / v1,
     * then set (u2u1u0) -= quot(v1v0)
     */
    t = ul;
    quot = t / v1;
    rem = t - quot * v1;
    /* Test on v0 */
    t = SHL(rem) | u0;
    if (B == quot || (quot * v0) > t) {
        quot--;
        rem += v1;
        t = SHL(rem) | u0;
        if (rem < B && (quot * v0) > t) {
            quot--;
        }
    }

    /*
     * multiply and subtract:
     * (u2u1u0)" = (u2u1u0)' - quot(v1v0)
     */
    uh = MSH(ul);           /* (0u2) */
    ul = SHL(ul) | u0;  /* (u1u0) */

    p0 = quot * v0;
    p1 = quot * v1;
    s = p0 + SHL(p1);
    ul -= s;
    borrow = (ul > ~s ? 1 : 0);
    uh -= MSH(p1) - borrow;
    if (0 != MSH(uh)) {
        /* add back v */
        quot--;
        ul += v;
        uh = 0;
    }

    /* quotient q1q0 */
    *q = SHL(q1) | LSH(quot);

    /* Remainder is in (u1u0) i.e. ul */
    *r = ul;

    /* quotient q2 (overflow) is returned */
    return q2;
}

#endif  /* HAVE_U64 */

/**
 * @brief Return greatest common divisor of two single limbs
 *
 * Returns gcd(x, y)
 *
 * Ref: Schneier 2nd ed, p245
 *
 * @param x single limb candidate #1
 * @param y single limb candidate #2
 *
 * @result return zero if x and y are zero, else gcd(x,y)
 */
limb_t sl_gcd(limb_t x, limb_t y) {
    limb_t g;

    if (x + y == 0) {
        return 0;    /* Error */
    }

    g = y;
    while (x > 0) {
        g = x;
        x = y % x;
        y = g;
    }
    return g;
}

/**
 * @brief Compute single limb exp = x^e mod m
 *
 * Computes exp = x^e mod m
 * Binary left-to-right method
 *
 * @param exp pointer to limb to receive result
 * @param x single limb x (base)
 * @param e single limb e (exponent)
 * @param m single limb m (modulus)
 *
 * @result zero on success (always!?)
 */
int sl_modexp(limb_t *exp, limb_t x, limb_t e, limb_t m) {
    limb_t mask;
    limb_t y;   /* Temp variable */

    /* Find most significant bit in e */
    for (mask = DIGMSB; mask > 0; mask >>= 1) {
        if (e & mask) {
            break;
        }
    }

    y = x;

    for (mask >>= 1; mask > 0; mask >>= 1) {
        sl_modmul(&y, y, y, m);     /* y = (y^2) % m */
        if (e & mask) {
            sl_modmul(&y, y, x, m);    /* y = (y*x) % m*/
        }
    }

    *exp = y;
    return 0;
}

/**
 * @brief Compute single limb inverse inv = u^(-1) % v
 *
 * Computes inv = u^(-1) % v
 * Ref: Knuth Algorithm X Vol 2 p 342
 * ignoring u2, v2, t2 and avoiding negative numbers
 *
 * @param inv pointer to limb to receive result
 * @param u single limb to inverse
 * @param v single limb modulus
 *
 * @result zero on success (always!?)
 */
int sl_modinv(limb_t *inv, limb_t u, limb_t v) {
    limb_t u1, u3, v1, v3, t1, t3, q, w;
    int iter = 1;

    /* Step X1. Initialize */
    u1 = 1;
    u3 = u;
    v1 = 0;
    v3 = v;

    /* Step X2. */
    while (v3 != 0) {
        /* Step X3. */
        q = u3 / v3;    /* Divide and */
        t3 = u3 % v3;
        w = q * v1; /* "Subtract" */
        t1 = u1 + w;
        /* Swap */
        u1 = v1;
        v1 = t1;
        u3 = v3;
        v3 = t3;
        iter = -iter;
    }

    if (iter < 0) {
        *inv = v - u1;
    } else {
        *inv = u1;
    }

    return 0;
}

/**
 * @brief Compute single limb a = (x * y) % mod
 *
 * Computes a = (x * y) % m
 *
 * @param a pointer to single limb to receive result
 * @param x single limb factor 1
 * @param y single limb factor 2
 * @param m single limb modulus
 *
 * @result zero on success (always!?)
 */
int sl_modmul(limb_t *a, limb_t x, limb_t y, limb_t m) {
    static IGRAPH_THREAD_LOCAL limb_t pp[2];

    /* pp[] = x * y */
    sl_mul(pp, x, y);

    /* *a = pp[] % m */
    *a = bn_mod_limb(pp, m, 2);

    /* Clean temp */
    pp[0] = pp[1] = 0;
    return 0;
}

#if HAVE_U64
/**
 * @brief Compute double limb product of two single limbs
 *
 * Computes p[] = x * y
 * where p is two limbs (double precision) and x, y are single
 * limbs. Use double precision natively supported on this machine.
 *
 * @param p pointer to an array of two limbs receiving the result
 * @param x single limb factor #1
 * @param y single limb factor #2
 *
 * @result zero on success (always)
 */
int sl_mul(limb_t p[2], limb_t x, limb_t y) {
    dlimb_t dd;

    dd = (dlimb_t)x * y;
    p[0] = (limb_t)dd;
    p[1] = (limb_t)(dd >> 32);
    return 0;
}

#else

/**
 * @brief Compute double limb product of two single limbs
 *
 * Computes p[] = x * y
 * Source: Arbitrary Precision Computation
 * http://numbers.computation.free.fr/Constants/constants.html
 *
 * The limbs x and y are split in halves and the four products
 * x1*y1, x0*y1, x1*y0 and x0*y0 are added shifting them to
 * their respective least significant bit position:
 * p[1] = x1*y1 + high(x0*y1 + x1*y0) + ch << 16 + cl
 * p[0] = x0*y0 + low(x0*y1 + x1*y0) << 16
 * ch = carry from adding x0*y1 + x1*y0
 * cl = carry from adding low(x0*y1 + x1*y0) << 16 to p[0]
 *
 * @param p pointer to an array of two limbs receiving the result
 * @param x single limb factor #1
 * @param y single limb factor #2
 *
 * @result zero on success (always)
 */
int sl_mul(limb_t p[2], limb_t x, limb_t y) {
    limb_t x0, y0, x1, y1;
    limb_t t, u, carry;

    /*
     * Split each x,y into two halves
     *   x = x0 + B*x1
     *   y = y0 + B*y1
     * where B = 2^16, half the limb size
     * Product is
     *   xy = x0y0 + B(x0y1 + x1y0) + B^2(x1y1)
     */
    x0 = LSH(x);
    x1 = MSH(x);
    y0 = LSH(y);
    y1 = MSH(y);

    /* Compute low part (w/o carry) */
    p[0] = x0 * y0;

    /* middle part */
    t = x0 * y1;
    u = x1 * y0;
    t += u;
    carry = (t < u ? 1 : 0);

    /*
     * The carry will go to high half of p[1],
     * and the high half of t will go into the
     * into low half of p[1]
     */
    carry = SHL(carry) + MSH(t);

    /* add low half of t to high half of p[0] */
    t = SHL(t);
    p[0] += t;
    if (p[0] < t) {
        carry++;
    }

    p[1] = x1 * y1 + carry;

    return 0;
}

#endif  /* HAVE_U64 */

/**
 * @brief Compute division of big number by a "half digit"
 *
 * Computes q[] = u[] / v, also returns r = u[] % v
 * where q, a are big number integers of nlimb limbs each,
 * and d, r are single limbs
 *
 * Using bit-by-bit method from MSB to LSB,
 * so v must be <= HALFMASK
 *
 * According to "Principles in PGP by Phil Zimmermann"
 *
 * @param q pointer to an array of limbs to receive the result
 * @param u pointer to an array of limbs (dividend)
 * @param v single limb (actually half limb) divisor
 * @param nlimb number of limbs in the arrays
 *
 * @result returns remainder of the division
 */
limb_t bn_div_hdig(limb_t q[], limb_t u[], limb_t v, count_t nlimb) {
    limb_t mask = DIGMSB;
    limb_t r = 0;
    if (v > HALFMASK) {
        igraph_errorf("bn_div_hdig called with v:%x", __FILE__,
                      __LINE__, (int) v);
    }

    if (0 == nlimb) {
        return 0;
    }
    if (0 == v) {
        return 0;    /* Divide by zero error */
    }

    /* Initialize quotient */
    bn_zero(q, nlimb);

    /* Work from MSB to LSB */
    while (nlimb > 0) {
        /* Multiply remainder by 2 */
        r <<= 1;

        /* Look at current bit */
        if (u[nlimb - 1] & mask) {
            r++;
        }
        if (r >= v) {
            /* Remainder became greater than divisor */
            r -= v;
            q[nlimb - 1] |= mask;
        }

        /* next bit */
        mask >>= 1;
        if (0 != mask) {
            continue;
        }

        /* next limb */
        --nlimb;
        mask = DIGMSB;
    }
    return r;
}

/**
 * @brief Compute single limb remainder of bignum % single limb
 *
 * Computes r = u[] % v
 * where a is a big number integer of nlimb
 * and r, v are single limbs, using bit-by-bit
 * method from MSB to LSB.
 *
 * Ref:
 *   Derived from principles in PGP by Phil Zimmermann
 * Note:
 *   This method will only work until r <<= 1 overflows.
 *   i.e. for d < DIGMSB, but we keep HALF_DIGIT
 *   limit for safety, and also because we don't
 *   have a 32nd bit.
 *
 * @param u pointer to big number to divide
 * @param v single limb (actually half limb) modulus
 * @param nlimb number of limbs in the array
 *
 * @result returns remainder of the division
 */
limb_t bn_mod_hdig(limb_t u[], limb_t v, count_t nlimb) {
    limb_t mask;
    limb_t r;

    if (0 == nlimb) {
        return 0;
    }
    if (0 == v) {
        return 0;    /* Divide by zero error */
    }

    if (v > HALFMASK) {
        igraph_errorf("bn_mod_hdig called with v:%x", __FILE__,
                      __LINE__, (int) v);
    }

    /* Work from left to right */
    mask = DIGMSB;
    r = 0;
    while (nlimb > 0) {
        /* Multiply remainder by 2 */
        r <<= 1;

        /* Look at current bit */
        if (u[nlimb - 1] & mask) {
            r++;
        }

        if (r >= v)
            /* Remainder became greater than divisor */
        {
            r -= v;
        }

        /* next bit */
        mask >>= 1;
        if (0 != mask) {
            continue;
        }

        /* next limb */
        --nlimb;
        mask = DIGMSB;
    }
    return r;
}

/**
 * @brief Addition of two bignum arrays
 *
 * Computes w[] = u[] + v[]
 * where w, u, v are big number integers of nlimb limbs each.
 * Returns carry, i.e. w[nlimb], as 0 or 1.
 *
 * Ref: Knuth Vol 2 Ch 4.3.1 p 266 Algorithm A.
 *
 * @param w pointer to array of limbs to receive the result
 * @param u pointer to array of limbs (addend #1)
 * @param v pointer to array of limbs (addend #2)
 * @param nlimb number of limbs in the arrays
 *
 * @result returns the carry, i.e. w[nlimb], as 0 or 1
 */
limb_t bn_add(limb_t w[], limb_t u[], limb_t v[], count_t nlimb) {
    limb_t carry;
    count_t j;

    for (j = 0, carry = 0; j < nlimb; j++) {
        /*
         * add limbs w[j] = u[j] + v[j] + carry;
         * set carry = 1 if carry (overflow) occurs
         */
        w[j] = u[j] + carry;
        carry = (w[j] < carry ? 1 : 0);

        w[j] = w[j] + v[j];
        if (w[j] < v[j]) {
            carry++;
        }
    }

    /* w[n] = carry */
    return carry;
}

/**
 * @brief Subtraction of two bignum arrays
 *
 * Calculates w[] = u[] - v[] where u[] >= v[]
 * w, u, v are big number integers of nlimb limbs each
 * Returns 0 if ok, or 1 if v was greater than u.
 *
 * Ref: Knuth Vol 2 Ch 4.3.1 p 267 Algorithm S.
 *
 * @param w pointer to array of limbs to receive the result
 * @param u pointer to array of limbs (minuend)
 * @param v pointer to array of limbs (subtrahend)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on success, 1 if v was greater than u
 */
limb_t bn_sub(limb_t w[], limb_t u[], limb_t v[], count_t nlimb) {
    limb_t borrow;
    count_t j;

    for (j = 0, borrow = 0; j < nlimb; j++) {
        /*
         * Subtract limbs w[j] = u[j] - v[j] - borrow;
         * set borrow = 1 if borrow occurs
         */
        w[j] = u[j] - borrow;
        borrow = (w[j] > ~borrow ? 1 : 0);

        w[j] = w[j] - v[j];
        if (w[j] > ~v[j]) {
            borrow++;
        }
    }

    /* borrow should be 0, if u >= v */
    return borrow;
}

/**
 * @brief Product of two bignum arrays
 *
 * Computes product w[] = u[] * v[]
 * where u, v are big number integers of nlimb each
 * and w is a big number integer of 2*nlimb limbs.
 *
 * Ref: Knuth Vol 2 Ch 4.3.1 p 268 Algorithm M.
 *
 * @param w pointer to array of limbs to receive the result
 * @param u pointer to array of limbs (factor #1)
 * @param v pointer to array of limbs (factor #2)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on success (always!?)
 */
int bn_mul(limb_t w[], limb_t u[], limb_t v[], count_t nlimb) {
    limb_t t[2];
    limb_t carry;
    count_t i, j, m, n;

    m = n = nlimb;

    /* zero result */
    bn_zero(w, 2 * nlimb);

    for (j = 0; j < n; j++) {
        /* zero multiplier? */
        if (0 == v[j]) {
            w[j + m] = 0;
            continue;
        }
        /* Initialize i */
        carry = 0;
        for (i = 0; i < m; i++) {
            /*
             * Multiply and add:
             * t = u[i] * v[j] + w[i+j] + carry
             */
            sl_mul(t, u[i], v[j]);

            t[0] += carry;
            if (t[0] < carry) {
                t[1]++;
            }
            t[0] += w[i + j];
            if (t[0] < w[i + j]) {
                t[1]++;
            }

            w[i + j] = t[0];
            carry = t[1];
        }
        w[j + m] = carry;
    }

    return 0;
}

/**
 * @brief Shift left a bignum by a number of bits (less than LIMBBITS)
 *
 * Computes a[] = b[] << x
 * Where a and b are big number integers of nlimb each.
 * The shift count must be less than LIMBBITS
 *
 * @param a pointer to array of limbs to receive the result
 * @param b pointer to array of limbs to shift left
 * @param x number of bits to shift (must be less than LIMBBITS)
 * @param nlimb number of limbs in the arrays
 *
 * @result returns a single limb "carry", i.e. bits that came out left
 */
limb_t bn_shl(limb_t a[], limb_t b[], count_t x, count_t nlimb) {
    count_t i, y;
    limb_t carry, temp;

    if (0 == nlimb) {
        return 0;
    }

    if (0 == x) {
        /* no shift at all */
        if (a != b) {
            bn_copy(a, b, nlimb);
        }
        return 0;
    }

    /* check shift amount */
    if (x >= LIMBBITS) {
        igraph_errorf("bn_shl() called with x >= %d", __FILE__,
                      __LINE__, LIMBBITS);
        return 0;
    }

    y = LIMBBITS - x;
    carry = 0;
    for (i = 0; i < nlimb; i++) {
        temp = b[i] >> y;
        a[i] = (b[i] << x) | carry;
        carry = temp;
    }

    return carry;
}

/**
 * @brief Shift right a bignum by a number of bits (less than LIMBBITS)
 *
 * Computes a[] = b[] >> x
 * Where a and b are big number integers of nlimb each.
 * The shift count must be less than LIMBBITS
 *
 * @param a pointer to array of limbs to receive the result
 * @param b pointer to array of limbs to shift right
 * @param x number of bits to shift (must be less than LIMBBITS)
 * @param nlimb number of limbs in the arrays
 *
 * @result returns a single limb "carry", i.e. bits that came out right
 */
limb_t bn_shr(limb_t a[], limb_t b[], count_t x, count_t nlimb) {
    count_t i, y;
    limb_t carry, temp;

    if (0 == nlimb) {
        return 0;
    }

    if (0 == x) {
        /* no shift at all */
        if (a != b) {
            bn_copy(a, b, nlimb);
        }
        return 0;
    }

    /* check shift amount */
    if (x >= LIMBBITS) {
        igraph_errorf("bn_shr() called with x >= %d", __FILE__,
                      __LINE__, LIMBBITS);
    }

    y = LIMBBITS - x;
    carry = 0;
    i = nlimb;
    while (i-- > 0) {
        temp = b[i] << y;
        a[i] = (b[i] >> x) | carry;
        carry = temp;
    }

    return carry;
}

/**
 * @brief Check a quotient for overflow
 *
 * Returns 1 if quot is too big,
 * i.e. if (quot * Vn-2) > (b.rem + Uj+n-2)
 * Returns 0 if ok
 *
 * @param quot quotient under test
 * @param rem remainder
 * @param
 *
 * @result zero on success
 */
static int quot_overflow(limb_t quot, limb_t rem, limb_t v, limb_t u) {
    limb_t t[2];

    sl_mul(t, quot, v);
    if (t[1] < rem) {
        return 0;
    }
    if (t[1] > rem) {
        return 1;
    }
    if (t[0] > u) {
        return 1;
    }

    return 0;
}

/**
 * @brief Compute quotient and remainder of bignum division
 *
 * Computes quotient q[] = u[] / v[]
 * and remainder r[] = u[] % v[]
 * where q, r, u are big number integers of ulimb limbs,
 * and the divisor v of vlimb limbs.
 *
 * Ref: Knuth Vol 2 Ch 4.3.1 p 272 Algorithm D.
 *
 * @param q pointer to array of limbs to receive quotient
 * @param r pointer to array of limbs to receive remainder
 * @param u pointer to array of limbs (dividend)
 * @param ulimb number of limbs in the q, r, u arrays
 * @param v pointer to array of limbs (divisor)
 * @param vlimb number of limbs in the v array
 *
 * @result zero on success, LIMBASK on division by zero
 */
int bn_div(limb_t q[], limb_t r[], limb_t u[], limb_t v[],
           count_t ulimb, count_t vlimb) {
    static IGRAPH_THREAD_LOCAL limb_t qq[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t uu[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t vv[BN_MAXSIZE];
    limb_t mask;
    limb_t overflow;
    limb_t quot;
    limb_t rem;
    limb_t t[2];
    limb_t *ww;
    count_t n, m, i, j, shift;
    int ok, cmp;

    /* find size of v */
    n = bn_sizeof(v, vlimb);

    /* Catch special cases */
    if (0 == n) {
        return (int) LIMBMASK;    /* Error: divide by zero */
    }

    if (1 == n) {
        /* Use short division instead */
        r[0] = bn_div_limb(q, u, v[0], ulimb);
        return 0;
    }

    /* find size of u */
    m = bn_sizeof(u, ulimb);

    if (m < n) {
        /* v > u: just set q = 0 and r = u */
        bn_zero(q, ulimb);
        bn_copy(r, u, ulimb);
        return 0;
    }

    if (m == n) {
        /* u and v are the same length: compare them */
        cmp = bn_cmp(u, v, (unsigned int)n);
        if (0 == cmp) {
            /* v == u: set q = 1 and r = 0 */
            bn_limb(q, 1, ulimb);
            bn_zero(r, ulimb);
            return 0;
        }
        if (cmp < 0) {
            /* v > u: set q = 0 and r = u */
            bn_zero(q, ulimb);
            bn_copy(r, u, ulimb);
            return 0;
        }
    }

    /* m greater than or equal to n */
    m -= n;

    /* clear quotient qq */
    bn_zero(qq, ulimb);

    /*
     * Normalize v: requires high bit of v[n-1] to be set,
     * so find most significant bit, then shift left
     */
    mask = DIGMSB;
    for (shift = 0; shift < LIMBBITS; shift++) {
        if (v[n - 1] & mask) {
            break;
        }
        mask >>= 1;
    }

    /* normalize vv from v */
    overflow = bn_shl(vv, v, shift, n);

    /* copy normalized dividend u into remainder uu */
    overflow = bn_shl(uu, u, shift, n + m);

    /* new limb u[m+n] */
    t[0] = overflow;

    j = m + 1;
    while (j-- > 0) {
        /* quot = (b * u[j+n] + u[j+n-1]) / v[n-1] */
        ok = 0;

        /* This is Uj+n */
        t[1] = t[0];
        t[0] = uu[j + n - 1];

        overflow = sl_div(&quot, &rem, t, vv[n - 1]);

        if (overflow) {
            /* quot = b */
            quot = LIMBMASK;
            rem = uu[j + n - 1] + vv[n - 1];
            if (rem < vv[n - 1]) {
                ok = 1;
            }
        }
        if (0 == ok && quot_overflow(quot, rem, vv[n - 2], uu[j + n - 2])) {
            /* quot * v[n-2] > b * rem + u[j+n-2] */
            quot--;
            rem += vv[n - 1];
            if (rem >= vv[n - 1])
                if (quot_overflow(quot, rem, vv[n - 2], uu[j + n - 2])) {
                    quot--;
                }
        }

        /* multiply and subtract vv[] * quot */
        ww = &uu[j];

        if (0 == quot) {
            overflow = 0;
        } else {
            /* quot is non zero */
            limb_t tt[2];
            limb_t borrow;

            for (i = 0, borrow = 0; i < n; i++) {
                sl_mul(tt, quot, vv[i]);
                ww[i] -= borrow;
                borrow = (ww[i] > ~borrow ? 1 : 0);

                ww[i] -= tt[0];
                if (ww[i] > ~tt[0]) {
                    borrow++;
                }
                borrow += tt[1];
            }

            /*
             * w[n] is not in array w[0..n-1]:
             * subtract final borrow
             */
            overflow = t[1] - borrow;
        }

        /* test for remainder */
        if (overflow) {
            quot--;
            /* add back if mul/sub was negative */
            overflow = bn_add(ww, ww, vv, n);
        }

        qq[j] = quot;

        /* u[j+n] for next round */
        t[0] = uu[j + n - 1];
    }

    /* clear uu[] limbs from n to n+m */
    for (j = n; j < m + n; j++) {
        uu[j] = 0;
    }

    /* denormalize remainder */
    bn_shr(r, uu, shift, n);

    /* copy quotient */
    bn_copy(q, qq, n + m);

    /* clear temps */
    bn_zero(qq, n);
    bn_zero(uu, n);
    bn_zero(vv, n);
    return 0;
}

/**
 * @brief Compute remainder of bignum division (modulo)
 *
 * Calculates r[] = u[] % v[]
 * where r, v are big number integers of length vlimb
 * and u is a big number integer of length ulimb.
 * r may overlap v.
 *
 * Note that r here is only vlimb long,
 * whereas in bn_div it is ulimb long.
 *
 * Use remainder from bn_div function.
 *
 * @param r pointer to array of limbs to receive remainder
 * @param u pointer to array of limbs (dividend)
 * @param ulimb number of limbs in the u array
 * @param v pointer to array of limbs (divisor)
 * @param vlimb number of limbs in the r and v array
 *
 * @result zero on success, LIMBASK on division by zero
 */
limb_t bn_mod(limb_t r[], limb_t u[], count_t ulimb, limb_t v[], count_t vlimb) {
    static IGRAPH_THREAD_LOCAL limb_t qq[2 * BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t rr[2 * BN_MAXSIZE];
    limb_t d0;

    /* rr[] = u[] % v[n] */
    d0 = (limb_t) bn_div(qq, rr, u, v, ulimb, vlimb);

    /* copy vlimb limbs of remainder */
    bn_copy(r, rr, vlimb);

    /* zero temps */
    bn_zero(rr, ulimb);
    bn_zero(qq, ulimb);

    return d0;
}

/**
 * @brief Compute greatest common divisor
 *
 * Computes g = gcd(x, y)
 * Reference: Schneier
 *
 * @param g pointer to array of limbs to receive the gcd
 * @param x pointer to array of limbs (candidate #1)
 * @param y pointer to array of limbs (candidate #2)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on succes (always)
 */
int bn_gcd(limb_t g[], limb_t x[], limb_t y[], count_t nlimb) {
    static IGRAPH_THREAD_LOCAL limb_t yy[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t xx[BN_MAXSIZE];

    bn_copy(xx, x, nlimb);
    bn_copy(yy, y, nlimb);

    /* g = y */
    bn_copy(g, yy, nlimb);

    /* while (x > 0) { */
    while (0 != bn_cmp_limb(xx, 0, nlimb)) {
        /* g = x */
        bn_copy(g, xx, nlimb);
        /* x = y % x */
        bn_mod(xx, yy, nlimb, xx, nlimb);
        /* y = g */
        bn_copy(yy, g, nlimb);
    }

    bn_zero(xx, nlimb);
    bn_zero(yy, nlimb);

    /* gcd is left in g */
    return 0;
}

/**
 * @brief Compute modular exponentiation of bignums
 *
 * Computes y[] = (x[]^e[]) % m[]
 * Binary MSB to LSB method
 *
 * @param y pointer to array of limbs to receive the result
 * @param x pointer to array of limbs (base)
 * @param e pointer to array of limbs (exponent)
 * @param m pointer to array of limbs (modulus)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on success, -1 on error (nlimb is zero)
 */
int bn_modexp(limb_t y[], limb_t x[], limb_t e[], limb_t m[], count_t nlimb) {
    limb_t mask;
    count_t n;

    if (nlimb == 0) {
        return -1;
    }

    /* Find second-most significant bit in e */
    n = bn_sizeof(e, nlimb);
    for (mask = DIGMSB; 0 != mask; mask >>= 1) {
        if (e[n - 1] & mask) {
            break;
        }
    }
    /* next bit, because we start off with y[] == x[] */
    mask >>= 1;
    if (0 == mask) {
        mask = DIGMSB;
        n--;
    }

    /* y[] = x[] */
    bn_copy(y, x, nlimb);

    while (n > 0) {
        /* y[] = (y[] ^ 2) % m[] */
        bn_modmul(y, y, y, m, nlimb);

        if (e[n - 1] & mask)
            /* y[] = (y[] * x[]) % m[] */
        {
            bn_modmul(y, y, x, m, nlimb);
        }

        /* next bit */
        mask >>= 1;
        if (0 == mask) {
            mask = DIGMSB;
            n--;
        }
    }

    return 0;
}

/**
 * @brief Compute modular product of two bignums
 *
 * Computes a[] = (x[] * y[]) % m[]
 * where a, x, y and m are big numbers of nlimb length
 *
 * @param a pointer to array of limbs to receive the result
 * @param x pointer to array of limbs (factor #1)
 * @param y pointer to array of limbs (factor #2)
 * @param m pointer to array of limbs (modulus)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on success, LIMBMASK if m was zero (division by zero)
 */
limb_t bn_modmul(limb_t a[], limb_t x[], limb_t y[], limb_t m[], count_t nlimb) {
    static IGRAPH_THREAD_LOCAL limb_t pp[2 * BN_MAXSIZE];
    limb_t d0;

    /* pp[] = x[] * y[] (NB: double size pp[]) */
    bn_mul(pp, x, y, nlimb);

    /* a[] = pp[] % m[] */
    d0 = bn_mod(a, pp, 2 * nlimb, m, nlimb);

    /* zero temp */
    bn_zero(pp, 2 * nlimb);

    return d0;
}

/**
 * @brief Compute modular inverse
 *
 * Computes inv[] = u[]^(-1) % v[]
 * Ref: Knuth Algorithm X Vol 2 p 342
 * ignoring u2, v2, t2 and avoiding negative numbers.
 *
 * @param inv pointer to array of limbs receiving the result
 * @param u pointer to array of limbs (candidate)
 * @param v pointer to array of limbs (modulus)
 * @param nlimb number of limbs in the arrays
 *
 * @result zero on success
 */
int bn_modinv(limb_t inv[], limb_t u[], limb_t v[], count_t nlimb) {
    /* Allocate temp variables */
    static IGRAPH_THREAD_LOCAL limb_t u1[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t u3[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t v1[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t v3[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t t1[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t t3[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t q[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t w[2 * BN_MAXSIZE];
    int iter;

    /* Step X1. Initialize */
    bn_limb(u1, 1, nlimb);  /* u1 = 1 */
    bn_limb(v1, 0, nlimb);  /* v1 = 0 */
    bn_copy(u3, u, nlimb);  /* u3 = u */
    bn_copy(v3, v, nlimb);  /* v3 = v */

    /* remember odd/even iterations */
    iter = 1;

    /* Step X2. Loop while v3 != 0 */
    while (0 != bn_cmp_limb(v3, 0, nlimb)) {
        /* Step X3. Divide and "Subtract" */
        /* q = u3 / v3, t3 = u3 % v3 */
        bn_div(q, t3, u3, v3, nlimb, nlimb);
        /* w = q * v1 */
        bn_mul(w, q, v1, nlimb);
        /* t1 = u1 + w */
        bn_add(t1, u1, w, nlimb);

        /* Swap u1 <= v1 <= t1 */
        bn_copy(u1, v1, nlimb);
        bn_copy(v1, t1, nlimb);

        /* Swap u3 <= v3 <= t3 */
        bn_copy(u3, v3, nlimb);
        bn_copy(v3, t3, nlimb);

        iter ^= 1;
    }

    if (iter) {
        bn_copy(inv, u1, nlimb);    /* inv = u1 */
    } else {
        bn_sub(inv, v, u1, nlimb);    /* inv = v - u1 */
    }

    /* clear temp vars */
    bn_zero(u1, nlimb);
    bn_zero(v1, nlimb);
    bn_zero(t1, nlimb);
    bn_zero(u3, nlimb);
    bn_zero(v3, nlimb);
    bn_zero(t3, nlimb);
    bn_zero(q, nlimb);
    bn_zero(w, 2 * nlimb);

    return 0;
}

/**
 * @brief Compute square root (and fraction) of a bignum
 *
 * Compute q[] = sqrt(u[]),
 * where q and u are big number integers of nlimb limbs
 *
 * Method according to sqrt.html of 2001-08-15:
 * Act on bytes from MSB to LSB, counting the number of times
 * that we can subtract consecutive odd numbers starting with
 * 1, 3, 5. Just uses add, subtract, shift and comparisons.
 *
 * The pointer r can be NULL if caller is not interested in
 * the (partial) fraction.
 *
 * @param q pointer to array of limbs to receive the result (integer)
 * @param r pointer to array of limbs to receive the result (fraction)
 * @param u pointer to array of limbs (square)
 * @param rlimb number of limbs in the q and r arrays
 * @param ulimb number of limbs in the u array
 *
 * @result zero on success
 */
int bn_sqrt(limb_t q[], limb_t r[], limb_t u[], count_t rlimb, count_t ulimb) {
    static IGRAPH_THREAD_LOCAL limb_t step[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t accu[BN_MAXSIZE];
    static IGRAPH_THREAD_LOCAL limb_t w[2 * BN_MAXSIZE];
    limb_t d;
    count_t m, n;
    count_t shift;

    bn_zero(q, ulimb);
    bn_limb(step, 1, BN_MAXSIZE);
    bn_limb(accu, 0, BN_MAXSIZE);
    n = bn_sizeof(u, ulimb);

    /* determine first non-zero byte from MSB to LSB */
    if (0 != (u[n - 1] >> 24)) {
        shift = 32;
    } else if (0 != (u[n - 1] >> 16)) {
        shift = 24;
    } else if (0 != (u[n - 1] >> 8)) {
        shift = 16;
    } else {
        shift = 8;
    }

    m = 1;
    while (n-- > 0) {
        while (shift > 0) {
            /* shift accu one byte left */
            bn_shl(accu, accu, 8, m + 1);

            /* shift for next byte from u[] */
            shift -= 8;
            accu[0] |= (u[n] >> shift) & 0xff;

            /* digit = 0 */
            d = 0;
            /* subtract consecutive odd numbers step[] until overflow */
            for (d = 0; bn_cmp(step, accu, m + 1) <= 0; d++) {
                bn_sub(accu, accu, step, m + 1);
                bn_add_limb(step, step, 2, m + 1);
            }

            /* put digit into result */
            bn_shl(q, q, 4, m);
            q[0] |= d;

            /* step[] = 2 * q[] * 16 + 1 */
            bn_shl(step, q, 5, m + 1);
            bn_add_limb(step, step, 1, m + 1);
        }
        shift = 32;
        if (0 == (n & 1)) {
            m++;
        }
    }

    /* Caller does not want to know the fraction? */
    if (NULL == r) {
        return 0;
    }

    /* nothing left to do if remainder is zero */
    if (0 == bn_cmp_limb(accu, 0, ulimb)) {
        bn_zero(r, rlimb);
        return 0;
    }

    /* Start off with the integer part */
    bn_zero(w, 2 * BN_MAXSIZE);
    bn_copy(w, q, ulimb);

    n = rlimb * (LIMBBITS / 4);
    while (n-- > 0) {
        /* shift accu one byte left */
        bn_shl(accu, accu, 8, rlimb);

        /* subtract consecutive odd numbers step[] until overflow */
        for (d = 0; bn_cmp(step, accu, rlimb) <= 0; d++) {
            bn_sub(accu, accu, step, rlimb);
            bn_add_limb(step, step, 2, rlimb);
        }

        /* put digit into result */
        bn_shl(w, w, 4, rlimb);
        w[0] |= d;

        /* step[] = 2 * w[] * 16 + 1 */
        bn_shl(step, w, 5, rlimb);
        bn_add_limb(step, step, 1, rlimb);
    }

    /* copy remainder */
    bn_copy(r, w, rlimb);
    return 0;
}
