/*****************************************************************************
 *  Entropy - Emerging Network To Reduce Orwellian Potency Yield
 *
 *  Copyright (C) 2005 Juergen Buchmueller <pullmoll@t-online.de>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
 *  $Id: bignum.h,v 1.6 2005/08/11 17:57:39 pullmoll Exp $
 *****************************************************************************/
#ifndef _bignum_h_
#define _bignum_h_

#include "config.h"
#ifdef HAVE_STDINT_H
    #include <stdint.h>
#else
    #ifdef HAVE_SYS_INT_TYPES_H
        #include <sys/int_types.h>
    #else
        #include "pstdint.h"
    #endif
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef NULL
    #define NULL 0
#endif

#ifndef O_BINARY
    #define O_BINARY 0
#endif

#ifndef HAVE_U64
    #define HAVE_U64 1
#endif

/* up to 512 limbs (512 * 32 = 16384 bits) numbers */
/* BN_MAXSIZE used to be 512 here, allowing us to go up to 512*32 = 16384 bits.
 * However, this has caused compilation problems with clang 7.3 (unless
 * compiling with -O2 -g). Since it is unlikely that we'll need that many bits,
 * I have changed this to 128, which still yields 4096 bits of precision but
 * does not cause problems with clang -- TN, 2016-04-18 */
#define BN_MAXSIZE 128
#define LIMBBITS 32
#define LIMBMASK 0xfffffffful
#define HALFMASK 0x0000fffful
#define DIGMSB 0x80000000ul
#define DIGLSB 0x00000001ul

typedef uint32_t count_t;
typedef uint16_t half_t;
typedef uint32_t limb_t;
#if HAVE_U64
    typedef uint64_t dlimb_t;
#endif

/* less significant half limb */
#define LSH(d)  ((half_t)(d))
/* more significant half limb */
#define MSH(d)  ((limb_t)(d)>>16)
/* shift left half limb */
#define SHL(d)  ((limb_t)(d)<<16)

/* single limb functions */
limb_t sl_div(limb_t *q, limb_t *r, limb_t u[2], limb_t v);
limb_t sl_gcd(limb_t x, limb_t y);
int sl_modexp(limb_t *exp, limb_t x, limb_t n, limb_t d);
int sl_modinv(limb_t *inv, limb_t u, limb_t v);
int sl_modmul(limb_t *a, limb_t x, limb_t y, limb_t m);
int sl_mul(limb_t p[2], limb_t x, limb_t y);

/* big number functions (max. MAXSIZE limbs) */
void bn_zero(limb_t a[], count_t nlimb);
void bn_limb(limb_t a[], limb_t d, count_t nlimb);
void bn_copy(limb_t a[], limb_t b[], count_t nlimb);
count_t bn_sizeof(limb_t a[], count_t nlimb);
int bn_cmp_limb(limb_t a[], limb_t b, count_t nlimb);
int bn_cmp(limb_t a[], limb_t b[], count_t nlimb);

/* big number to hex, decimal, binary */
const char *bn2x(limb_t a[], count_t nlimb);
const char *bn2d(limb_t a[], count_t nlimb);
const char *bn2f(limb_t a[], count_t alimb, limb_t b[], count_t blimb);
const char *bn2b(limb_t a[], count_t nlimb);

/* big number with single limb operations */
limb_t bn_add_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb);
limb_t bn_sub_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb);
limb_t bn_div_limb(limb_t q[], limb_t u[], limb_t v, count_t nlimb);
limb_t bn_mod_limb(limb_t u[], limb_t d, count_t nlimb);
limb_t bn_mul_limb(limb_t w[], limb_t u[], limb_t v, count_t nlimb);

/* big number with single limb <= HALFMASK operations */
limb_t bn_div_half(limb_t q[], limb_t u[], limb_t v, count_t nlimb);
limb_t bn_mod_half(limb_t a[], limb_t d, count_t nlimb);

/* big number operations */
limb_t bn_add(limb_t w[], limb_t u[], limb_t v[], count_t nlimb);
limb_t bn_sub(limb_t w[], limb_t u[], limb_t v[], count_t nlimb);
limb_t bn_shl(limb_t a[], limb_t b[], count_t x, count_t nlimb);
limb_t bn_shr(limb_t a[], limb_t b[], count_t x, count_t nlimb);
int bn_mul(limb_t w[], limb_t u[], limb_t v[], count_t nlimb);
int bn_div(limb_t q[], limb_t r[], limb_t u[], limb_t v[], count_t ulimb, count_t vlimb);
limb_t bn_mod(limb_t r[], limb_t u[], count_t ulimb, limb_t v[], count_t vlimb);
int bn_gcd(limb_t g[], limb_t x[], limb_t y[], count_t nlimb);
int bn_sqrt(limb_t g[], limb_t x[], limb_t y[], count_t rlimb, count_t nlimb);
int bn_modexp(limb_t y[], limb_t x[], limb_t e[], limb_t m[], count_t nlimb);
int bn_modinv(limb_t inv[], limb_t u[], limb_t v[], count_t nlimb);
limb_t bn_modmul(limb_t a[], limb_t x[], limb_t y[], limb_t m[], count_t nlimb);

#endif  /* !defined(_bignum_h_) */
