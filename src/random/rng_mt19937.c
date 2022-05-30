/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_random.h"

#include "igraph_memory.h"
#include "igraph_types.h"

#include "config.h" /* IGRAPH_THREAD_LOCAL */

#include <string.h> /* memset() */


#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;

typedef struct {
    unsigned long mt[N];
    int mti;
} igraph_i_rng_mt19937_state_t;

static igraph_uint_t igraph_rng_mt19937_get(void *vstate) {
    igraph_i_rng_mt19937_state_t *state = vstate;

    unsigned long k;
    unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

    if (state->mti >= N) {
        /* generate N words at one time */
        int kk;

        for (kk = 0; kk < N - M; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
        for (; kk < N - 1; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

        {
            unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
        }

        state->mti = 0;
    }

#undef MAGIC

    /* Tempering */

    k = mt[state->mti];
    k ^= (k >> 11);
    k ^= (k << 7) & 0x9d2c5680UL;
    k ^= (k << 15) & 0xefc60000UL;
    k ^= (k >> 18);

    state->mti++;

    return k;
}

static igraph_error_t igraph_rng_mt19937_seed(void *vstate, igraph_uint_t seed) {
    igraph_i_rng_mt19937_state_t *state = vstate;
    int i;

    memset(state, 0, sizeof(igraph_i_rng_mt19937_state_t));

    if (seed == 0) {
        seed = 4357;   /* the default seed is 4357 */
    }
    state->mt[0] = seed & 0xffffffffUL;

    for (i = 1; i < N; i++) {
        /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
           Ed. p.106 for multiplier. */
        state->mt[i] =
            (1812433253UL * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) +
             (unsigned long) i);
        state->mt[i] &= 0xffffffffUL;
    }

    state->mti = i;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_mt19937_init(void **state) {
    igraph_i_rng_mt19937_state_t *st;

    st = IGRAPH_CALLOC(1, igraph_i_rng_mt19937_state_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize MT19937 RNG");
    (*state) = st;

    igraph_rng_mt19937_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_mt19937_destroy(void *vstate) {
    igraph_i_rng_mt19937_state_t *state =
        (igraph_i_rng_mt19937_state_t*) vstate;
    IGRAPH_FREE(state);
}

/**
 * \var igraph_rngtype_mt19937
 * \brief The MT19937 random number generator.
 *
 * The MT19937 generator of Makoto Matsumoto and Takuji Nishimura is a
 * variant of the twisted generalized feedback shift-register
 * algorithm, and is known as the “Mersenne Twister” generator. It has
 * a Mersenne prime period of 2^19937 - 1 (about 10^6000) and is
 * equi-distributed in 623 dimensions. It has passed the diehard
 * statistical tests. It uses 624 words of state per generator and is
 * comparable in speed to the other generators. The original generator
 * used a default seed of 4357 and choosing \c s equal to zero in
 * \c igraph_rng_mt19937_seed() reproduces this. Later versions switched to
 * 5489 as the default seed, you can choose this explicitly via
 * \ref igraph_rng_seed() instead if you require it.
 *
 * </para><para>
 * For more information see,
 * Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A
 * 623-dimensionally equidistributed uniform pseudorandom number
 * generator”. ACM Transactions on Modeling and Computer Simulation,
 * Vol. 8, No. 1 (Jan. 1998), Pages 3–30
 *
 * </para><para>
 * The generator \c igraph_rngtype_mt19937 uses the second revision of the
 * seeding procedure published by the two authors above in 2002. The
 * original seeding procedures could cause spurious artifacts for some
 * seed values.
 *
 * </para><para>
 * This generator was ported from the GNU Scientific Library.
 */

const igraph_rng_type_t igraph_rngtype_mt19937 = {
    /* name= */      "MT19937",
    /* bits=  */     32,
    /* init= */      igraph_rng_mt19937_init,
    /* destroy= */   igraph_rng_mt19937_destroy,
    /* seed= */      igraph_rng_mt19937_seed,
    /* get= */       igraph_rng_mt19937_get,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0
};

/***** Default RNG, used upon igraph startup *****/

#define addr(a) (&a)

static igraph_i_rng_mt19937_state_t igraph_i_rng_default_state = {
    /* .mt = */ {
        /* Pre-calculated state for the default seed. It is used to prevent
         * igraph_rng_get_integer() from an infinite loop when the RNG is
         * used without seeding. */
        0x00001105, 0x9c4563fa, 0x25a828da, 0xef2ac805, 0x41162062, 0x2bc4c214,
        0x4f0945ea, 0x328058be, 0x3152b0fe, 0xb549c23f, 0x0ab7471b, 0xa51980b2,
        0x8dbdf57c, 0x72c448c3, 0x27ae8698, 0xdc0a7207, 0x353f1fa4, 0xd7523fc5,
        0x3dfe1f30, 0x8e42fe03, 0x5863bf79, 0x1d19c26d, 0x91310a17, 0xa9513760,
        0x1f5c4bc2, 0xae42b5a3, 0x56ebd19f, 0xe1894171, 0x514ad416, 0x9ea6fc30,
        0xf23641d8, 0xc9a32e86, 0x42e98799, 0x34dad719, 0x65223dff, 0x13a16359,
        0x6045d341, 0x20599865, 0x7caa2cff, 0xa429ae5d, 0xeb5ba2a3, 0xd78ec949,
        0x2e14045c, 0xd2bef477, 0xf55485f0, 0x33ace40c, 0x359464ea, 0x07f30a81,
        0xe9092e15, 0x64e2f4df, 0x9ca769c8, 0x822fd6e5, 0x4b326857, 0xc7133023,
        0xa0331cd6, 0xc863d3db, 0x66a92c70, 0xec9e01ce, 0x71e86b1b, 0x37b92b7d,
        0x00bf0d8d, 0x647bcdde, 0x9b1d9039, 0xda767a86, 0x0df283b9, 0xef3af93e,
        0x7e50fa53, 0x9c26a49d, 0xddb109ff, 0xcb14ccb1, 0xf19e0480, 0x4756e2f6,
        0x0d7ebabb, 0x240abf10, 0xf1ecf19a, 0x66f532a8, 0x85786df9, 0x013bb754,
        0xe4f74872, 0x526b0de4, 0xcbef08a9, 0x55976563, 0xe7aa71fc, 0x4f3970ee,
        0xe9a2759f, 0x254ce2e1, 0xfce6ec1b, 0xf01bfdcf, 0xb6714dd4, 0x06353bc7,
        0x866f14dd, 0x22179356, 0x6080274a, 0x5ea3a3f4, 0xc100cd07, 0xde2406f3,
        0x5c7d2d10, 0xe1f7e116, 0xd1d30aab, 0xe38f1cab, 0x066436ac, 0x357f9e41,
        0x26d1390b, 0x7b5f64be, 0x0ec5f6c3, 0xc87db658, 0xc0aca551, 0x85d41bc5,
        0x05fd74ef, 0x071c09b8, 0xb34b4e06, 0xdb91ec03, 0x2bde1c70, 0x32e928a1,
        0xe52230f7, 0xf95be4b7, 0xbf918f78, 0xde8be597, 0x6219c7da, 0x5b1d0cde,
        0xec726b73, 0xb4b353a9, 0xb03385f1, 0x44a7e45a, 0x91eecb63, 0xd9b426c2,
        0xc81793a3, 0x18adde9f, 0xf314ec3b, 0x2f332a99, 0xb72aafdf, 0x26ffa7b4,
        0xb9887e88, 0xec57c6f7, 0x98c612ca, 0x29aa716f, 0x94fc2853, 0xff46417e,
        0xd62ebbdb, 0xc5deb4c3, 0x30d8104c, 0x680b1a89, 0x174c4036, 0x88ec3bdd,
        0x2856f68b, 0xd60ca868, 0xf80cb5c9, 0xdad2d345, 0x7727d132, 0x8e0ed4b4,
        0x94a95264, 0x8cb018d5, 0x20a5dc6b, 0xcf5239d0, 0x6527bbd9, 0x6f1bb4d3,
        0x3370b976, 0x78f3522b, 0xe61de530, 0x91d8b8be, 0xe2817ecc, 0x759ccf4c,
        0x78e8ff03, 0x3572ae6d, 0x369226a5, 0x68da8cbe, 0xc0c4bf01, 0x68e56e71,
        0x6fae82d8, 0x93cec146, 0xeedba47e, 0x0319cafc, 0xe3b1f218, 0x1281f854,
        0x1577edd2, 0xbf543689, 0xaf36e887, 0xfabbea2a, 0x2f9253df, 0x72b66eae,
        0xbdfe52bf, 0x0ec8ca46, 0x59614454, 0x0b447340, 0x3a6eb8f8, 0xfd73b291,
        0x63339654, 0xb720cd44, 0xf2b5735a, 0x8fc023da, 0xa4e4bcf6, 0x92092103,
        0xd34a8f25, 0x350fd0bf, 0xaf2c931d, 0x112aa2fe, 0x33fe3cfa, 0xb4c5d967,
        0xc128d29f, 0xa0149453, 0x65b4ddbd, 0x4b2917f5, 0x3cb4080e, 0x6bb8ac51,
        0xca42cc5c, 0x36527948, 0x42696236, 0x7ca42f82, 0x9ac8da7f, 0x43961922,
        0x71a0a6a1, 0xdbed5df3, 0x9d7b8084, 0x452b6bb3, 0x90a1c010, 0x781869f1,
        0x97e33c88, 0x67d2bd4b, 0x2277490c, 0x7e7b3e97, 0xd937f80a, 0x3b10ad6a,
        0xb04825b0, 0x48812219, 0xc4d94c58, 0x6e14d3c8, 0x2a0d202f, 0x52afda6e,
        0x884f95af, 0x723aa326, 0x5c833e49, 0xd0141b4f, 0xa99e71e4, 0xdbc106a7,
        0x5e31639e, 0x652565a6, 0x26b87acf, 0xb02a3b98, 0xb3c4eeb0, 0xef4d6f29,
        0x9c305682, 0xdadca171, 0xed8fb4ec, 0x680f4a3e, 0xd87d02cf, 0x73654771,
        0xee402026, 0xcd837c90, 0xde8fd1f7, 0x73c56a3d, 0x8a6106a6, 0xe35163af,
        0x87275fd8, 0xfbca7bff, 0x6720c76a, 0xd293f136, 0x34e687e9, 0x410e50ee,
        0x5e7ed64d, 0x14cd38ff, 0x2b70f49f, 0x5ad29ac0, 0x312f582b, 0x4b09ccfe,
        0x088b58a3, 0x21de3458, 0xd515bfc2, 0x5f79f130, 0x31866261, 0x1151ba52,
        0x42326568, 0x2806347c, 0xdbed11fc, 0xb52291ac, 0x5a5b98b8, 0xe3ed4310,
        0xbf05a293, 0xff2fbd4a, 0x3320bfe3, 0xb5d130a6, 0x330bf5cc, 0x52d62695,
        0xc72e6d7e, 0x2854186c, 0x25336fb8, 0xe21c8cb5, 0x1d8aeaec, 0xfbfafc3b,
        0x40877b38, 0x60f41f9e, 0xd69591dd, 0x59245bb9, 0xc275a8bc, 0xe0f2cb80,
        0x834866d5, 0x0e7aa2fa, 0x447017ca, 0x47820740, 0x6af5a6cf, 0x20d10e71,
        0xe44c2cc1, 0xe74f7bb7, 0x1d7b2332, 0x95c8a5e9, 0x266139e7, 0x35a97854,
        0x978d6e56, 0xef237c57, 0x2cd50258, 0xdeb305ed, 0x3c4fb61c, 0x3eabd643,
        0xcf4a64a7, 0x88b179ed, 0x11cf0385, 0x043191b4, 0x1a73d140, 0x236ccf7d,
        0x9d62c28f, 0xbfef37e0, 0x9d76ff6a, 0xaa796d49, 0xde7242d9, 0xb4d10b45,
        0x78737347, 0xfd26f1e3, 0xd4ee4ea6, 0x4c9e5560, 0x37cf998d, 0xccef0aea,
        0x5f850037, 0xa60ffc99, 0x5eba9d73, 0xa4002147, 0xcabe0e87, 0xc85c5f63,
        0x971c0230, 0xc297a10b, 0xb237d17a, 0x8966ddab, 0x9dd3e601, 0x38c65b84,
        0xbdd7c06a, 0x7eee925f, 0x69060e6e, 0xe6261a24, 0x5c1831bd, 0xf44d3c87,
        0x0f778570, 0xf6a0968d, 0x43cf6564, 0x1fd90f38, 0xc44ffa78, 0xfdf6a6e8,
        0xe6109f19, 0xfb69b0a5, 0x9cc588e2, 0xce4ee1c5, 0x505c0a84, 0x5b9454e0,
        0xc118e72d, 0x3cccd48f, 0xbe7964d5, 0x92bcd93e, 0x3566d218, 0xf6a7bce5,
        0x468e9e2c, 0x78257e30, 0x96aa03c5, 0x7b89fdf4, 0xe40c501b, 0xe06072eb,
        0x27e37efc, 0x4499f6e1, 0xb5fc47d6, 0x35c4cc1b, 0x809bfb1f, 0x50b798ea,
        0x6bb21931, 0xbd0ea16b, 0x2105e0e9, 0x9f0d6e6a, 0xd53a3886, 0xefda7af8,
        0xcddfd987, 0xa7547695, 0x29ea9a15, 0x289708cc, 0x24dfa600, 0xc6127f85,
        0x10350764, 0x299c6ffb, 0x1682828f, 0x9f450af4, 0xddd4fa98, 0xbd60d3b2,
        0x924cb5fc, 0xe79abcc3, 0x584f394e, 0x2215e4ba, 0xe320c8f2, 0x520f41a6,
        0xc9b84775, 0xf02c5921, 0x45205dfe, 0xbe0b8e30, 0x2505dd50, 0x8bf02227,
        0x43064731, 0xffe2c789, 0xfe07950c, 0xaf2bd786, 0x273dacb0, 0x758f530d,
        0x6a51315a, 0x22ef2d86, 0x8e63ad7e, 0x938ecf8d, 0xa5546c0d, 0xb08baa8e,
        0x1f3336e0, 0xb7b28805, 0xb06c6b69, 0xb42fa5de, 0xef932d74, 0xbaa6a09c,
        0xdc03ee00, 0xae0183da, 0x310f9de4, 0x37e450a1, 0x5ca0fa33, 0xa0c57969,
        0x91d02be7, 0xbde1e00a, 0x4b39acda, 0xa93f671a, 0xa5d7862c, 0xd8188fdb,
        0x759359ee, 0x87cc6402, 0xb12775b8, 0xb2a7fe1b, 0xe1f1a297, 0x56625a1f,
        0xd5d99d92, 0x4f24c4f2, 0xcd8cc09d, 0x57ee8e15, 0xadb2c3a4, 0x93c5083f,
        0xb3cfe6d3, 0xaf42eb38, 0xdbded9a6, 0x90e82cde, 0x759f7092, 0xd32b16c6,
        0x1a926a81, 0x2c4a0fae, 0xf15e4f70, 0x62dee52a, 0xb2c16ec3, 0x7ddcfcf2,
        0x5a32d8ad, 0x48b789ab, 0x64bb4be2, 0xd0b66d60, 0xba3824e1, 0x4a190a62,
        0xaf2515e3, 0x777e0c9a, 0x55b6edfd, 0x886bc243, 0x393f6e7d, 0x49927e2a,
        0x037dcbd1, 0x716b4450, 0x31184ed1, 0xb47bf352, 0xf14cd06e, 0xfdd391e0,
        0x1ec00b6f, 0x0bebebac, 0xf6e708be, 0x05fe9974, 0xbfbda0a8, 0x3d715ef7,
        0xb64da859, 0x7a3920ce, 0x04c1ba93, 0xab4c48e8, 0x1e84003c, 0xffd83597,
        0x8c055950, 0x24272147, 0x113021f1, 0x23be5f04, 0x4410a284, 0x4f2c4d6a,
        0x33d5d029, 0xdce01320, 0x589a49c3, 0x62a5ed7f, 0x5e0122ac, 0xa0bc4538,
        0x87e95bda, 0x241fd631, 0xd186bc4f, 0xfb00f7f7, 0x43bf6940, 0x769f51a2,
        0x0f00724d, 0xc4734f60, 0x26af4f0f, 0x91e239ec, 0x9cc33ae8, 0xa7f37a55,
        0x9be8d557, 0xce59a98e, 0x83f75ba7, 0x8e257720, 0xbb773472, 0x7027a239,
        0x84fafa22, 0x1dbfd0ab, 0x3c05d883, 0xfec188bc, 0x65c22c69, 0x26372f17,
        0x9999e523, 0x253d1116, 0x1ed685c0, 0xa07686d3, 0x98a00b89, 0xd21ef2ec,
        0xeec0c161, 0xfb37bfc1, 0x08eb7ba2, 0x437f7b03, 0x0f2d9be4, 0xd8a7870f,
        0xe2b1b5d8, 0xe76bf484, 0xb724ba61, 0xc5f1862e, 0x005306e1, 0xab9721e6,
        0xb0ff6516, 0xb5669707, 0xb56e441d, 0x87ce7960, 0xbd1957d0, 0x29bd0a01,
        0x53f67d8d, 0x88437665, 0x8dcbd7cd, 0x4292edd6, 0x1a20e6ff, 0xcf929bc8,
        0x40c41c45, 0xf65b8d03, 0x0b93a330, 0x49e41421, 0x42a112d2, 0x04645a72,
        0x011cb32e, 0x1a78515b, 0x727ece1d, 0x7b184f43, 0x95cf9942, 0x64aab879,
        0xbf5b0192, 0xfdadb00b, 0x38f2bd64, 0xbcdb3eb1, 0x11f289dd, 0xe182ab70,
        0x1a73319f, 0x8459acfc, 0x17e73078, 0x1a6e599b, 0x52ba4f6b, 0x3be71117,
        0x75f00f59, 0x8b5427ff, 0xd5852e19, 0xbded1c8b, 0x66909557, 0x344ff339,
        0x1f4378c9, 0x8ac13a9a, 0xe0bf7846, 0x88ca6288, 0x305abcc2, 0xa41e4cdb,
        0x9d0374ef, 0xfa6af8d4, 0x8d403e27, 0x6a9d53ee, 0x618606a1, 0xe1cc3f77,
        0x0fb61f1c, 0x53444465, 0x714481ce, 0x3a280006, 0xe1f53aba, 0x093c2e5a,
        0x160275e0, 0xb92d63bf, 0xcf7380f1, 0x473163db, 0x37e21164, 0xb13b62d7,
        0x4424fd6d, 0x2c2aca01, 0x45033dcb, 0xdede7d19, 0x627047aa, 0x0653cbe0,
        0x6a4a51ca, 0x95a4ea82, 0xb50f06ec, 0x12251c53, 0xcf119a2d, 0xbcb67495
    },
    /* .mti = */ N
};

IGRAPH_THREAD_LOCAL igraph_rng_t igraph_i_rng_default = {
    addr(igraph_rngtype_mt19937),
    addr(igraph_i_rng_default_state),
    /* is_seeded = */ 1
};

#undef addr
#undef N
#undef M
