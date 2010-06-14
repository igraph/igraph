/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph_random.h"
#include "igraph_error.h"
#include "config.h"

#include <math.h>
#include <limits.h>
#include "igraph_math.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_memory.h"

/** 
 * \section about_rngs 
 * 
 * <section>
 * <title>About random numbers in igraph, use cases</title>
 * 
 * <para> 
 * Some algorithms in igraph, e.g. the generation of random graphs,
 * require random number generators (RNGs). Prior to version 0.6
 * igraph did not have a sophisticated way to deal with random number
 * generators at the C level, but this has changed. From version 0.6
 * different and multiple random number generators are supported.
 * </para>
 * </section>
 * 
 */

/** 
 * \section rng_use_cases
 * 
 * <section><title>Use cases</title>
 * 
 * <section><title>Normal (default) use</title>
 * <para> 
 * If the user does use any of the RNG functions explicitly, but calls
 * some of the randomized igraph functions, then a default RNG is set
 * up, the first time an igraph function needs random numbers. The
 * seed of this RNG is the output of the <code>time(0)</code> function
 * call, using the <code>time</code> function from the standard C
 * library. This ensures that igraph creates a different random graph,
 * easch time the C program is called.
 * </para>
 * 
 * <para> 
 * The created default generator is stored in the \ref
 * igraph_rng_default variable.
 * </para>
 * </section>
 * 
 * <section><title>Reproducible simulations</title>
 * <para> 
 * If reproducible results are needed, then the user should set the
 * seed of the default random number generator explixitly, using the 
 * \ref igraph_rng_seed() function on the default generator, \ref
 * igraph_rng_default. When setting the seed to the same number,
 * igraph generates exactly the same random graph (or series of random
 * graphs).
 * </para>
 * </section>
 * 
 * <section><title>Changing the default generator</title>
 * <para> 
 * TODO
 * </para>
 * </section>
 *
 * <section><title>Use multiple generators</title>
 * <para> 
 * TODO
 * </para>
 * </section>
 *
 * </section>
 */

/* ------------------------------------ */

typedef struct {
  int i, j;
  long int x[31];
} igraph_i_rng_glibc2_state_t;

unsigned long int igraph_i_rng_glibc2_get(int *i, int *j, int n, 
					  long int *x) {
  long int k;

  x[*i] += x[*j];
  k = (x[*i] >> 1) & 0x7FFFFFFF;
  
  (*i)++;
  if (*i == n) {
    *i = 0;
  }
  
  (*j)++ ;
  if (*j == n) {
    *j = 0;
  }

  return k;
}

unsigned long int igraph_rng_glibc2_get(void *vstate) {
  igraph_i_rng_glibc2_state_t *state = 
    (igraph_i_rng_glibc2_state_t*) vstate;
  return igraph_i_rng_glibc2_get(&state->i, &state->j, 31, state->x);
}

igraph_real_t igraph_rng_glibc2_get_real(void *state) {
  return igraph_rng_glibc2_get(state) / 2147483648.0;
}

/* this function is independent of the bit size */

void igraph_i_rng_glibc2_init(long int *x, int n, 
			      unsigned long int s) {
  int i;
  
  if (s==0) { s=1; }
  
  x[0] = s;
  for (i=1 ; i<n ; i++) {
    const long int h = s / 127773;
    const long int t = 16807 * (s - h * 127773) - h * 2836;
    if (t < 0) {
      s = t + 2147483647 ;
    } else { 
      s = t ;
    }
    
    x[i] = s ;
  }
}

int igraph_rng_glibc2_seed(void *vstate, unsigned long int seed) {
  igraph_i_rng_glibc2_state_t *state = 
    (igraph_i_rng_glibc2_state_t*) vstate;
  int i;
  
  igraph_i_rng_glibc2_init(state->x, 31, seed);
  
  state->i=3;
  state->j=0;
  
  for (i=0;i<10*31; i++) {
    igraph_rng_glibc2_get(state);
  }
  
  return 0;
}

int igraph_rng_glibc2_init(void **state) {
  igraph_i_rng_glibc2_state_t *st;

  st=igraph_Calloc(1, igraph_i_rng_glibc2_state_t);
  if (!st) {
    IGRAPH_ERROR("Cannot initialize RNG", IGRAPH_ENOMEM);
  }
  (*state)=st;

  igraph_rng_glibc2_seed(st, 0);
  
  return 0;
}

void igraph_rng_glibc2_destroy(void *vstate) {
  igraph_i_rng_glibc2_state_t *state = 
    (igraph_i_rng_glibc2_state_t*) vstate;
  igraph_Free(state);
}

/**
 * \var igraph_rngtype_glibc2
 * \brief The random number generator introduced in GNU libc 2
 * 
 * 
 */

igraph_rng_type_t igraph_rngtype_glibc2 = {
  /* name= */      "LIBC",
  /* min=  */      0,
  /* max=  */      RAND_MAX,
  /* init= */      igraph_rng_glibc2_init,
  /* destroy= */   igraph_rng_glibc2_destroy,
  /* seed= */      igraph_rng_glibc2_seed,
  /* get= */       igraph_rng_glibc2_get,
  /* get_real= */  igraph_rng_glibc2_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

/* ------------------------------------ */

typedef struct {
  unsigned long int x;
} igraph_i_rng_rand_state_t;

unsigned long int igraph_rng_rand_get(void *vstate) {
  igraph_i_rng_rand_state_t *state = vstate;
  state->x = (1103515245 * state->x + 12345) & 0x7fffffffUL;
  return state->x;
}

igraph_real_t igraph_rng_rand_get_real(void *vstate) {
  return igraph_rng_rand_get (vstate) / 2147483648.0 ;
}

int igraph_rng_rand_seed(void *vstate, unsigned long int seed) {
  igraph_i_rng_rand_state_t *state = vstate;
  state->x = seed;
  return 0;
}

int igraph_rng_rand_init(void **state) {
  igraph_i_rng_rand_state_t *st;

  st=igraph_Calloc(1, igraph_i_rng_rand_state_t);
  if (!st) {
    IGRAPH_ERROR("Cannot initialize RNG", IGRAPH_ENOMEM);
  }
  (*state)=st;
  
  igraph_rng_rand_seed(st, 0);
  
  return 0;
}

void igraph_rng_rand_destroy(void *vstate) {
  igraph_i_rng_rand_state_t *state = 
    (igraph_i_rng_rand_state_t*) vstate;
  igraph_Free(state);  
}  

/**
 * \var igraph_rngtype_rand
 * \brief The old BSD rand/stand random number generator
 */

igraph_rng_type_t igraph_rngtype_rand = {
  /* name= */      "RAND",
  /* min=  */      0,
  /* max=  */      0x7fffffffUL,
  /* init= */      igraph_rng_rand_init,
  /* destroy= */   igraph_rng_rand_destroy,
  /* seed= */      igraph_rng_rand_seed,
  /* get= */       igraph_rng_rand_get,
  /* get_real= */  igraph_rng_rand_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

/* ------------------------------------ */

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

unsigned long int igraph_rng_mt19937_get(void *vstate) {
  igraph_i_rng_mt19937_state_t *state = vstate;

  unsigned long k ;
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

igraph_real_t igraph_rng_mt19937_get_real(void *vstate) {
  return igraph_rng_mt19937_get (vstate) / 4294967296.0 ;
}

int igraph_rng_mt19937_seed(void *vstate, unsigned long int seed) {
  igraph_i_rng_mt19937_state_t *state = vstate;
  int i;

  if (seed == 0) {
    seed = 4357;   /* the default seed is 4357 */
  }
  state->mt[0]= seed & 0xffffffffUL;

  for (i = 1; i < N; i++) {
    /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
       Ed. p.106 for multiplier. */
    state->mt[i] =
      (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
    state->mt[i] &= 0xffffffffUL;
  }
  
  state->mti = i;
  return 0;
}

int igraph_rng_mt19937_init(void **state) {
  igraph_i_rng_mt19937_state_t *st;

  st=igraph_Calloc(1, igraph_i_rng_mt19937_state_t);
  if (!st) {
    IGRAPH_ERROR("Cannot initialize RNG", IGRAPH_ENOMEM);
  }
  (*state)=st;
  
  igraph_rng_mt19937_seed(st, 0);
  
  return 0;
}

void igraph_rng_mt19937_destroy(void *vstate) {
  igraph_i_rng_mt19937_state_t *state = 
    (igraph_i_rng_mt19937_state_t*) vstate;
  igraph_Free(state);  
}  

/** 
 * \var igraph_rngtype_mt19937
 * \brief The MT19937 random number generator
 */

igraph_rng_type_t igraph_rngtype_mt19937 = {
  /* name= */      "MT19937",
  /* min=  */      0,
  /* max=  */      0xffffffffUL,
  /* init= */      igraph_rng_mt19937_init,
  /* destroy= */   igraph_rng_mt19937_destroy,
  /* seed= */      igraph_rng_mt19937_seed,
  /* get= */       igraph_rng_mt19937_get,
  /* get_real= */  igraph_rng_mt19937_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

/* int main() { */
/*   igraph_i_rng_mt19937_state_t *state; */
/*   int i; */
/*   igraph_rng_mt19937_init((void*) (&state)); */
/*   printf("%i\n", state->mti); */
/*   for (i=0;i<N;i++) { */
/*     printf("%li ", state->mt[i]); */
/*   } */
/*   printf("\n"); */
/*   return 0; */
/* } */

igraph_i_rng_mt19937_state_t igraph_i_rng_default_state = {
  { 
    4357, -1673174022, 631777498, -282408955, 1091969122, 734315028,
    1326007786, 847272126, 827502846, -1253457345, 179783451,
    -1525055310, -1916930692, 1925466307, 665749144, -603295225,
    893329316, -682475579, 1040064304, -1908212221, 1482932089,
    488227437, -1859057129, -1454295200, 526142402, -1371359837,
    1458295199, -511098511, 1363858454, -1633223632, -231325224,
    -912052602, 1122600857, 886757145, 1696742911, 329343833,
    1615188801, 542742629, 2091527423, -1540772259, -346316125,
    -678508215, 773063772, -759237513, -179010064, 866968588, 898917610,
    133368449, -385274347, 1692595423, -1666750008, -2110794011,
    1261594711, -955043805, -1607263018, -932981797, 1722362992,
    -325189170, 1911057179, 934882173, 12520845, 1685835230,
    -1692561351, -629769594, 233997241, -281347778, 2119236179,
    -1675189091, -575600129, -887829327, -241302400, 1196876534,
    226409147, 604684048, -236129894, 1727345320, -2055705095, 20690772,
    -453556110, 1382747620, -873527127, 1435985251, -408260100,
    1329164526, -375229025, 625795809, -51975141, -266601009,
    -1234088492, 104152007, -2039540515, 571970390, 1619011402,
    1587782644, -1056912121, -568064269, 1551707408, -503848682,
    -774698325, -477160277, 107230892, 897556033, 651245835, 2069849278,
    247854787, -931285416, -1062427311, -2049696827, 100496623,
    119278008, -1286910458, -611193853, 735976560, 854141089,
    -450744073, -111418185, -1080979592, -561257065, 1645856730,
    1528630494, -328045709, -1263316055, -1338800655, 1151853658,
    -1846621341, -642505022, -937978973, 414047903, -216732613,
    791882393, -1221939233, 654288820, -1182237048, -329791753,
    -1731849526, 699036015, -1795413933, -12172930, -701580325,
    -975260477, 819466316, 1745558153, 390873142, -1997784099,
    676787851, -703813528, -133384759, -623717563, 1999098162,
    -1911630668, -1800842652, -1934616363, 547740779, -816694832,
    1697102809, 1864086739, 863025526, 2029212203, -434248400,
    -1848067906, -494829876, 1973210956, 2028535555, 896708205,
    915547813, 1759153342, -1060847871, 1759866481, 1873707736,
    -1815166650, -287595394, 52022012, -474877416, 310507604, 360181202,
    -1085000055, -1355356025, -88348118, 798118879, 1924558510,
    -1107406145, 248040006, 1499546708, 189035328, 980334840, -42749295,
    1664325204, -1222587068, -222989478, -1883233318, -1528513290,
    -1844895485, -750088411, 890228927, -1356033251, 288006910,
    872299770, -1262102169, -1054289249, -1609264045, 1706352061,
    1260984309, 1018431502, 1807264849, -901591972, 911374664,
    1114202678, 2091134850, -1698112897, 1133910306, 1906353825,
    -605200909, -1652850556, 1160473523, -1868447728, 2014865905,
    -1746715512, 1741864267, 578242828, 2122006167, -650643446,
    990948714, -1337449040, 1216422425, -992392104, 1846858696,
    705503279, 1387256430, -2008050257, 1916445478, 1552105033,
    -803988657, -1449233948, -608106841, 1580295070, 1696949670,
    649624271, -1339409512, -1278939472, -280137943, -1674553726,
    -623074959, -309349140, 1745832510, -662895921, 1936017265,
    -297787354, -847020912, -560999945, 1942317629, -1973352794,
    -481205329, -2027462696, -70616065, 1730201450, -762056394,
    887523305, 1091457262, 1585370701, 348993791, 728822943, 1523751616,
    825186347, 1258933502, 143349923, 568210520, -719994942, 1601827120,
    830890593, 290568786, 1110599016, 671495292, -605220356,
    -1256025684, 1515952312, -470990064, -1090149741, -13648566,
    857784291, -1244581722, 856421836, 1389766293, -953258626,
    676599916, 624127928, -501445451, 495643372, -67437509, 1082620728,
    1626611614, -694840867, 1495555001, -1032476484, -520959104,
    -2092407083, 242918138, 1148196810, 1199703872, 1794483919,
    550571633, -464769855, -414221385, 494609202, -1782012439,
    643906023, 900298836, -1752338858, -282887081, 752157272,
    -558692883, 1011856924, 1051448899, -817208153, -2001634835,
    298779525, 70357428, 443797824, 594333565, -1654472049, -1074841632,
    -1653145750, -1434882743, -562937127, -1261368507, 2020832071,
    -47779357, -722579802, 1285444960, 936352141, -856749334,
    1602551863, -1508901735, 1589288307, -1543495353, -893514105,
    -933470365, -1759772112, -1030250229, -1304964742, -1989747285,
    -1647057407, 952523652, -1109933974, 2129564255, 1762004590,
    -433710556, 1545089469, -196264825, 259491184, -157247859,
    1137665380, 534318904, -1001391496, -34167064, -435118311,
    -76959579, -1664775966, -833691195, 1348209284, 1536447712,
    -1055332563, 1020056719, -1099340587, -1833117378, 895930904,
    -156779291, 1183751724, 2015723056, -1767242811, 2072641012,
    -468955109, -530550037, 669220604, 1150940897, -1241757738,
    902089755, -2137261281, 1354209514, 1806833969, -1123114645,
    554033385, -1626509718, -717604730, -270894344, -840967801,
    -1487636843, 703240725, 680986828, 618636800, -971866235, 271910756,
    698118139, 377651855, -1622865164, -573244776, -1117727822,
    -1840466436, -409289533, 1481587022, 571860154, -484390670,
    1376731558, -910669963, -265529055, 1159749118, -1106538960,
    621141328, -1947196889, 1124484913, -1914999, -33057524,
    -1356081274, 658353328, 1972327181, 1783705946, 586100102,
    -1906070146, -1819357299, -1521193971, -1333024114, 523450080,
    -1213036539, -1335071895, -1271945762, -275567244, -1163485028,
    -603722240, -1375632422, 823107044, 937709729, 1554053683,
    -1597671063, -1848628249, -1109270518, 1262071002, -1455462630,
    -1512602068, -669478949, 1972591086, -2016648190, -1322814024,
    -1297613285, -504257897, 1449286175, -707158638, 1327809778,
    -846413667, 1475251733, -1380793436, -1815803841, -1278220589,
    -1354568904, -606152282, -1863832354, 1973383314, -752150842,
    445803137, 743051182, -245477520, 1658774826, -1295946045,
    2111634674, 1513281709, 1219987883, 1689996258, -793350816,
    -1170725663, 1243155042, -1356524061, 2004749466, 1438051837,
    -2006203837, 960458365, 1234337322, 58575825, 1902855248, 823676625,
    -1266945198, -246624146, -36466208, 515902319, 200010668,
    -152631106, 100571508, -1078091608, 1030840055, -1236424615,
    2050564302, 79805075, -1421063960, 511967292, -2607721, -1945806512,
    606544199, 288367089, 599678724, 1141940868, 1328303466, 869650473,
    -589294816, 1486506435, 1655041407, 1577132716, -1598274248,
    -2014749734, 606066225, -779699121, -83822601, 1136617792,
    1990152610, 251687501, -999076000, 649023247, -1847445012,
    -1664927000, -1477215659, -1679239849, -832984690, -2080941145,
    -1910147296, -1149815694, 1881645625, -2063926750, 499110059,
    1007016067, -20870980, 1707224169, 641150743, -1717967581,
    624759062, 517375424, -1602844973, -1734341751, -769723668,
    -289357471, -80232511, 149650338, 1132428035, 254647268, -660109553,
    -491670056, -412355452, -1222329759, -974027218, 5441249,
    -1416158746, -1325439722, -1251567865, -1251064803, -2016511648,
    -1122412592, 700254721, 1408662925, -2008844699, -1916020787,
    1116925398, 438363903, -812475448, 1086594117, -161772285,
    194224944, 1239684129, 1117852370, 73685618, 18658094, 444092763,
    1920912925, 2065190723, -1781556926, 1688909945, -1084554862,
    -38948853, 955432292, -1126482255, 301107677, -511530128, 443756959,
    -2074497796, 401027192, 443439515, 1387941739, 1004998935,
    1978666841, -1957419009, -712692199, -1108534133, 1720751447,
    877654841, 524515529, -1967048038, -524322746, -2000002424,
    811252930, -1541518117, -1660717841, -93652780, -1925169625,
    1788695534, 1636173473, -506708105, 263593756, 1396982885,
    1900315086, 975699974, -504022342, 154938970, 369260000,
    -1188207681, -814513935, 1194419163, 937562468, -1321508137,
    1143274861, 741001729, 1157840331, -555844327, 1651525546,
    106154976, 1783255498, -1784354174, -1257306388, 304421971,
    -820930003, -1128893291     
  },
  N
};

#undef N
#undef M

/* ------------------------------------ */

#ifndef USING_R

/** 
 * \var igraph_rng_default
 * The default igraph random number generator
 */

igraph_rng_t igraph_rng_default = { 
  &igraph_rngtype_mt19937,
  &igraph_i_rng_default_state,
  /* def= */ 1
};

#endif

/* ------------------------------------ */

#ifdef USING_R

double  unif_rand(void);
double  norm_rand(void);
double  Rf_rgeom(double);
double  Rf_rbinom(double, double);

int igraph_rng_R_init(void **state) {
  IGRAPH_ERROR("R RNG error, unsupported function called",
	       IGRAPH_EINTERNAL);
  return 0;
}

void igraph_rng_R_destroy(void *state) {
  igraph_error("R RNG error, unsupported function called",
	       __FILE__, __LINE__, IGRAPH_EINTERNAL);
}

int igraph_rng_R_seed(void *state) {
  IGRAPH_ERROR("R RNG error, unsupported function called",
	       IGRAPH_EINTERNAL);
  return 0;
}

unsigned long int igraph_rng_R_get(void *state) {
  return unif_rand() * 0x7FFFFFFFUL;
}

igraph_real_t igraph_rng_R_get_real(void *state) {
  return unif_rand();
}

igraph_real_t igraph_rng_R_get_norm(void *state) {
  return norm_rand();
}

igraph_real_t igraph_rng_R_get_geom(void *state, igraph_real_t p) {
  return Rf_rgeom(p);
}
 
igraph_real_t igraph_rng_R_get_binom(void *state, long int n,
				     igraph_real_t p) {
  return Rf_rbinom(n, p);
}

igraph_rng_type_t igraph_rngtype_R = {
  /* name= */      "GNU R",
  /* min=  */      0,
  /* max=  */      0x7FFFFFFFUL,
  /* init= */      igraph_rng_R_init,
  /* destroy= */   igraph_rng_R_destroy,
  /* seed= */      igraph_rng_R_seed,
  /* get= */       igraph_rng_R_get,
  /* get_real= */  igraph_rng_R_get_real,
  /* get_norm= */  igraph_rng_R_get_norm,
  /* get_geom= */  igraph_rng_R_get_geom,
  /* get_binom= */ igraph_rng_R_get_binom
};

#endif

/* ------------------------------------ */

double igraph_norm_rand(igraph_rng_t *rng);
double igraph_rgeom(igraph_rng_t *rng, double p);
double igraph_rbinom(igraph_rng_t *rng, double nin, double pp);

/** 
 * \function igraph_rng_init
 * Initialize a random number generator
 */

int igraph_rng_init(igraph_rng_t *rng, const igraph_rng_type_t *type) {
  rng->type=type;
  IGRAPH_CHECK(rng->type->init(&rng->state));
  return 0;
}

/** 
 * \function igraph_rng_destroy
 * Deallocate memory associated with a random number generator
 */

void igraph_rng_destroy(igraph_rng_t *rng) {
  rng->type->destroy(rng->state);
}

/**
 * \function igraph_rng_seed
 * Set the seed of a random number generator
 */
int igraph_rng_seed(igraph_rng_t *rng, unsigned long int seed) {
  const igraph_rng_type_t *type=rng->type;
  rng->def=0;
  IGRAPH_CHECK(type->seed(rng->state, seed));
  return 0;
}

/** 
 * \function igraph_rng_max 
 * Query the maximum possible integer for a random number generator
 */

unsigned long int igraph_rng_max(igraph_rng_t *rng) {
  const igraph_rng_type_t *type=rng->type;
  return type->max;
}

/**
 * \function igraph_rng_min
 * Query the minimum possible integer for a random number generator
 */

unsigned long int igraph_rng_min(igraph_rng_t *rng) {
  const igraph_rng_type_t *type=rng->type;
  return type->min;
}

/** 
 * \function igraph_rng_name
 * Query the type of a random number generator
 */

const char *igraph_rng_name(igraph_rng_t *rng) {
  const igraph_rng_type_t *type=rng->type;
  return type->name;
}

/** 
 * \function igraph_rng_get_integer
 * Generate an integer random number from an interval
 */

long int igraph_rng_get_integer(igraph_rng_t *rng,
				long int l, long int h) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_real) {
    return (long int)(type->get_real(rng->state)*(h-l+1)+l);
  } else if (type->get) {
    unsigned long int max=type->max;
    return (long int)(type->get(rng->state))/
      ((double)max+1)*(h-l+1)+l;
  }
  IGRAPH_ERROR("Internal random generator error", IGRAPH_EINTERNAL);
  return 0;
}

/** 
 * \function igraph_rng_get_normal
 * Normally distributed random numbers
 */

igraph_real_t igraph_rng_get_normal(igraph_rng_t *rng, 
				    igraph_real_t m, igraph_real_t s) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_norm) {
    return type->get_norm(rng->state)*s+m;
  } else {
    return igraph_norm_rand(rng)*s+m;
  }
}

/** 
 * \function igraph_rng_get_unif
 * Generate real, uniform random numbers from an interval
 */

igraph_real_t igraph_rng_get_unif(igraph_rng_t *rng, 
				  igraph_real_t l, igraph_real_t h) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_real) {
    return type->get_real(rng->state)*(h-l)+l;
  } else if (type->get) {
    unsigned long int max=type->max;
    return type->get(rng->state)/((double)max+1)*(double)(h-l)+l;
  }
  IGRAPH_ERROR("Internal random generator error", IGRAPH_EINTERNAL);
  return 0;  
}

/** 
 * \function igraph_rng_get_unif01
 * Generate real, uniform random number from the unit interval
 */

igraph_real_t igraph_rng_get_unif01(igraph_rng_t *rng) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_real) {
    return type->get_real(rng->state);
  } else if (type->get) {
    unsigned long int max=type->max;
    return type->get(rng->state)/((double)max+1);
  }
  IGRAPH_ERROR("Internal random generator error", IGRAPH_EINTERNAL);
  return 0;  
}

/** 
 * \function igraph_rng_get_geom
 * Generate geometrically distributed random numbers
 */

igraph_real_t igraph_rng_get_geom(igraph_rng_t *rng, igraph_real_t p) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_geom) {
    return type->get_geom(rng->state, p);
  } else {
    return igraph_rgeom(rng, p);
  }
}

/** 
 * \function igraph_rng_get_binom
 * Generate binomially distributed random numbers
 */

igraph_real_t igraph_rng_get_binom(igraph_rng_t *rng, long int n, 
				   igraph_real_t p) {
  const igraph_rng_type_t *type=rng->type;
  if (type->get_binom) {
    return type->get_binom(rng->state, n, p);
  } else {
    return igraph_rbinom(rng, n, p);
  }
}

unsigned long int igraph_rng_get_int31(igraph_rng_t *rng) {
  const igraph_rng_type_t *type=rng->type;
  unsigned long int max=type->max;
  if (type->get && max==0x7FFFFFFFUL) {
    return type->get(rng->state);
  } else if (type->get_real) {
    return type->get_real(rng->state)*0x7FFFFFFFUL;
  } else { 
    return igraph_rng_get_unif01(rng)*0x7FFFFFFFUL;
  }
}

int igraph_rng_inited = 0;

#ifndef HAVE_EXPM1
#ifndef USING_R			/* R provides a replacement */
/* expm1 replacement */
static double expm1 (double x)
{
    if (fabs(x) < M_LN2)
    {
        /* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */

        double i = 1.0;
        double sum = x;
        double term = x / 1.0;

        do
        {
            term *= x / ++i;
            sum += term;
        }
        while (fabs(term) > fabs(sum) * 2.22e-16);
      
        return sum;
    }

    return expl(x) - 1.0L;
}
#endif
#endif

#ifndef HAVE_RINT
#ifndef USING_R			/* R provides a replacement */
/* rint replacement */
static double rint (double x)
{
   return ( (x<0.) ? -floor(-x+.5) : floor(x+.5) );
}
#endif
#endif

#ifndef HAVE_RINTF
static float rintf (float x)
{
   return ( (x<(float)0.) ? -(float)floor(-x+.5) : (float)floor(x+.5) );
}
#endif

/*
 * \ingroup internal
 * 
 * This function appends the rest of the needed random number to the 
 * result vector.
 */

int igraph_random_sample_alga(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			      igraph_integer_t length) {
  igraph_real_t N=h-l+1;
  igraph_real_t n=length;
  
  igraph_real_t top=N-n;
  igraph_real_t Nreal=N;
  igraph_real_t S=0;
  igraph_real_t V, quot;
  
  l=l-1;

  while (n>=2) {
    V=RNG_UNIF01();
    S=1;
    quot=top/Nreal;
    while (quot>V) {
      S+=1;
      top=-1.0+top;
      Nreal=-1.0+Nreal;
      quot=(quot*top)/Nreal;
    }
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
    Nreal=-1.0+Nreal; n=-1+n;
  }
  
  S=floor(round(Nreal)*RNG_UNIF01());
  l+=S+1;
  igraph_vector_push_back(res, l);	/* allocated */
  
  return 0;
}

/**
 * \ingroup nongraph
 * \function igraph_random_sample
 * \brief Generates an increasing random sequence of integers.
 * 
 * </para><para>
 * This function generates an incresing sequence of random integer
 * numbers from a given interval. The algorithm is taken literally
 * from Jeffrey Scott Vitter: 'An Efficient Algorithm for Sequential
 * Random Sampling', ACM Transactions on Mathematical Software, 13/1,
 * 58--67. This method can be used for generating numbers from a
 * \em very large interval, it is primilarly created for randomly
 * selecting some edges from the sometimes huge set of possible edges
 * in a large graph.
 * \param res Pointer to an initialized vector, this will hold the
 *        result. It will be resized to the proper size.
 * \param l The lower limit of the generation interval (inclusive).
 * \param h The upper limit of the generation interval (inclusive).
 * \param length The number of random integers to generate.
 * \return Error code.
 *
 * Time complexity: according to the referenced paper, the expected
 * running time is O(length).
 */

int igraph_random_sample(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			 igraph_integer_t length) {
  igraph_real_t N=h-l+1;
  igraph_real_t n=length;
  int retval;

  igraph_real_t nreal=length;
  igraph_real_t ninv=1.0/nreal;
  igraph_real_t Nreal=N;
  igraph_real_t Vprime;
  igraph_real_t qu1=-n+1+N;
  igraph_real_t qu1real=-nreal+1.0+Nreal;
  igraph_real_t negalphainv=-13;
  igraph_real_t threshold=-negalphainv*n;
  igraph_real_t S;
  
  igraph_vector_clear(res);
  IGRAPH_CHECK(igraph_vector_reserve(res, length));  

  RNG_BEGIN();
  
  Vprime=exp(log(RNG_UNIF01())*ninv);

  while (n>1 && threshold < N) {
    igraph_real_t X, U;
    igraph_real_t limit, t;
    igraph_real_t negSreal, y1, y2, top, bottom;
    igraph_real_t nmin1inv=1.0/(-1.0+nreal);
    while (1) {
      while(1) {
	X=Nreal*(-Vprime+1.0);
	S=floor(X);
	if (S==0) { S=1; }
	if (S <qu1) { break; }
	Vprime = exp(log(RNG_UNIF01())*ninv);
      }
      U=RNG_UNIF01();
      negSreal=-S;
      
      y1=exp(log(U*Nreal/qu1real)*nmin1inv);
      Vprime=y1*(-X/Nreal+1.0)*(qu1real/(negSreal+qu1real));
      if (Vprime <= 1.0) { break; }
      
      y2=1.0;
      top=-1.0+Nreal;
      if (-1+n > S) {
	bottom=-nreal+Nreal; 
	limit=-S+N;
      } else {
	bottom=-1.0+negSreal+Nreal;
	limit=qu1;
      }
      for (t=-1+N; t>=limit; t--) {
	y2=(y2*top)/bottom;
	top=-1.0+top;
	bottom=-1.0+bottom;
      }
      if (Nreal/(-X+Nreal) >= y1*exp(log(y2)*nmin1inv)) {
	Vprime=exp(log(RNG_UNIF01())*nmin1inv);
	break;
      }
      Vprime=exp(log(RNG_UNIF01())*ninv);
    }
        
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
    N=-S+(-1+N);   Nreal=negSreal+(-1.0+Nreal);
    n=-1+n;   nreal=-1.0+nreal; ninv=nmin1inv;
    qu1=-S+qu1; qu1real=negSreal+qu1real;
    threshold=threshold+negalphainv;
  }
  
  if (n>1) {
    retval=igraph_random_sample_alga(res, l, h, n);
  } else {
    retval=0;
    S=floor(N*Vprime);
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
  }

  RNG_END();
  
  return retval;
}
  
#ifdef USING_R

#else

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *	double qnorm5(double p, double mu, double sigma,
 *		      int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the quantile function for the normal distribution.
 *
 *	For small to moderate probabilities, algorithm referenced
 *	below is used to obtain an initial approximation which is
 *	polished with a final Newton step.
 *
 *	For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution,
 *	Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2004  The R Development Core Team
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef _MSC_VER
#  define ML_POSINF IGRAPH_INFINITY
#  define ML_NEGINF -IGRAPH_INFINITY
#  define ML_NAN    IGRAPH_NAN
#else
#  define ML_POSINF	(1.0 / 0.0)
#  define ML_NEGINF	((-1.0) / 0.0)
#  define ML_NAN		(0.0 / 0.0)
#endif

#define ML_ERROR(x)	/* nothing */
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define i1mach		Rf_i1mach
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr

	/* Chebyshev Series */

int	chebyshev_init(double*, int, double);
double	chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

void	gammalims(double*, double*);
double	lgammacor(double); /* log(gamma) correction */
double  stirlerr(double);  /* Stirling expansion "error" */

double	lfastchoose(double, double);

double  bd0(double, double);

/* Consider adding these two to the API (Rmath.h): */
double	dbinom_raw(double, double, double, double, int);
double	dpois_raw (double, double, int);
double  pnchisq_raw(double, double, double, double, double, int);

int	i1mach(int);

/* From toms708.c */
void bratio(double a, double b, double x, double y, 
	    double *w, double *w1, int *ierr);


#endif /* MATHLIB_PRIVATE_H */


	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? (p) : (1 - (p)))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (1 - (p)) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (1 - (p)))/* [log](1-p) */

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

/*till 1.8.x:
 * #define R_DT_val(x)	R_D_val(R_D_Lval(x))
 * #define R_DT_Cval(x)	R_D_val(R_D_Cval(x)) */
#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
/* ==   R_DT_log when we already "know" log_p == TRUE :*/

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }

double igraph_qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    if (p == R_DT_0)	return ML_NEGINF;
    if (p == R_DT_1)	return ML_POSINF;
    R_Q_P01_check(p);

    if(sigma  < 0)	ML_ERR_return_NAN;
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

	/* r = min(p, 1-p) < 0.075 */
	if (q > 0)
	    r = R_DT_CIv(p);/* 1-p */
	else
	    r = p_;/* = R_DT_Iv(p) ^=  p */

	r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

double fsign(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
        return x + y;
#endif
    return ((y >= 0) ? fabs(x) : -fabs(x));
}

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}

int imin2(int x, int y)
{
    return (x < y) ? x : y;
}

#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#elif HAVE_WORKING_FINITE
/* include header needed to define finite() */
#  ifdef HAVE_IEEE754_H
#   include <ieee754.h>         /* newer Linuxen */
#  else
#   ifdef HAVE_IEEEFP_H
#    include <ieeefp.h>         /* others [Solaris], .. */
#   endif
#  endif
# define R_FINITE(x)    finite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif

int R_finite(double x)
{
#ifdef HAVE_WORKING_ISFINITE
    return isfinite(x);
#elif HAVE_WORKING_FINITE
    return finite(x);
#else
/* neither finite nor isfinite work. Do we really need the AIX exception? */
# ifdef _AIX
#  include <fp.h>
     return FINITE(x);
# elif defined(_MSC_VER)
     return _finite(x);
#else
    return (!isnan(x) & (x != 1/0.0) & (x != -1.0/0.0));
# endif
#endif
}

int R_isnancpp(double x)
{
   return (isnan(x)!=0);
}

#ifdef __cplusplus
  int R_isnancpp(double); /* in arithmetic.c */
#  define ISNAN(x)     R_isnancpp(x)
#else
#  define ISNAN(x)     (isnan(x)!=0)
#endif

double igraph_norm_rand(igraph_rng_t *rng) {
  
  double u1;

#define BIG 134217728 /* 2^27 */
  /* unif_rand() alone is not of high enough precision */
  u1 = igraph_rng_get_unif01(rng);
  u1 = (int)(BIG*u1) + igraph_rng_get_unif01(rng);
  return igraph_qnorm5(u1/BIG, 0.0, 1.0, 1, 0);
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2002 the R Development Core Team
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double exp_rand(void);
 *
 *  DESCRIPTION
 *
 *    Random variates from the standard exponential distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1972).
 *    Computer methods for sampling from the exponential and
 *    normal distributions.
 *    Comm. ACM, 15, 873-882.
 */

double igraph_exp_rand(igraph_rng_t *rng)
{
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 8) is determined by q[n-1] = 1.0 */
    /* within standard precision */
    const double q[] =
    {
	0.6931471805599453,
	0.9333736875190459,
	0.9888777961838675,
	0.9984959252914960,
	0.9998292811061389,
	0.9999833164100727,
	0.9999985691438767,
	0.9999998906925558,
	0.9999999924734159,
	0.9999999995283275,
	0.9999999999728814,
	0.9999999999985598,
	0.9999999999999289,
	0.9999999999999968,
	0.9999999999999999,
	1.0000000000000000
    };
    double a, u, ustar, umin;
    int i;

    a = 0.;
    /* precaution if u = 0 is ever returned */
    u = igraph_rng_get_unif01(rng);
    while(u <= 0.0 || u >= 1.0) u = igraph_rng_get_unif01(rng);
    for (;;) {
	u += u;
	if (u > 1.0)
	    break;
	a += q[0];
    }
    u -= 1.;

    if (u <= q[0])
	return a + u;

    i = 0;
    ustar = igraph_rng_get_unif01(rng);
    umin = ustar;
    do {
        ustar = igraph_rng_get_unif01(rng);
	if (ustar < umin)
	    umin = ustar;
	i++;
    } while (u > q[i]);
    return a + umin * q[0];
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Development Core Team
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rpois(double lambda)
 *
 *  DESCRIPTION
 *
 *    Random variates from the Poisson distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1982).
 *    Computer generation of Poisson deviates
 *    from modified normal distributions.
 *    ACM Trans. Math. Software 8, 163-179.
 */

#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

#define repeat for(;;)

#define FALSE 0
#define TRUE  1
#define M_1_SQRT_2PI    0.398942280401432677939946059934     /* 1/sqrt(2pi) */

double igraph_rpois(igraph_rng_t *rng, double mu)
{
    /* Factorial Table (0:9)! */
    const double fact[10] =
    {
	1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static int l, m;

    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = FALSE;

    if (!R_FINITE(mu))
	ML_ERR_return_NAN;

    if (mu <= 0.)
	return 0.;

    big_mu = mu >= 10.;
    if(big_mu)
	new_big_mu = FALSE;

    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

	if (big_mu) {
	    new_big_mu = TRUE;
	    /* Case A. (recalculation of s,d,l	because mu has changed):
	     * The poisson probabilities pk exceed the discrete normal
	     * probabilities fk whenever k >= m(mu).
	     */
	    muprev = mu;
	    s = sqrt(mu);
	    d = 6. * mu * mu;
	    big_l = floor(mu - 1.1484);
	    /* = an upper bound to m(mu) for all mu >= 10.*/
	}
	else { /* Small mu ( < 10) -- not using normal approx. */

	    /* Case B. (start new table and calculate p0 if necessary) */

	    /*muprev = 0.;-* such that next time, mu != muprev ..*/
	    if (mu != muprev) {
		muprev = mu;
		m = imax2(1, (int) mu);
		l = 0; /* pp[] is already ok up to pp[l] */
		q = p0 = p = exp(-mu);
	    }

	    repeat {
		/* Step U. uniform sample for inversion method */
	        u = igraph_rng_get_unif01(rng);
		if (u <= p0)
		    return 0.;

		/* Step T. table comparison until the end pp[l] of the
		   pp-table of cumulative poisson probabilities
		   (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
		if (l != 0) {
		    for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
			if (u <= pp[k])
			    return (double)k;
		    if (l == 35) /* u > pp[35] */
			continue;
		}
		/* Step C. creation of new poisson
		   probabilities p[l..] and their cumulatives q =: pp[k] */
		l++;
		for (k = l; k <= 35; k++) {
		    p *= mu / k;
		    q += p;
		    pp[k] = q;
		    if (u <= q) {
			l = k;
			return (double)k;
		    }
		}
		l = 35;
	    } /* end(repeat) */
	}/* mu < 10 */

    } /* end {initialize persistent vars} */

/* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * igraph_norm_rand(rng);/* norm_rand() ~ N(0,1), standard normal */

    if (g >= 0.) {
	pois = floor(g);
	/* Step I. immediate acceptance if pois is large enough */
	if (pois >= big_l)
	    return pois;
	/* Step S. squeeze acceptance */
	fk = pois;
	difmuk = mu - fk;
	u = igraph_rng_get_unif01(rng); /* ~ U(0,1) - sample */
	if (d * u >= difmuk * difmuk * difmuk)
	    return pois;
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
	   because one might have exited in step I or S
	   */
        muprev2 = mu;
	omega = M_1_SQRT_2PI / s;
	/* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
	 * approximations to the discrete normal probabilities fk. */

	b1 = one_24 / mu;
	b2 = 0.3 * b1 * b1;
	c3 = one_7 * b1 * b2;
	c2 = b2 - 15. * c3;
	c1 = b1 - 6. * b2 + 45. * c3;
	c0 = 1. - b1 + 3. * b2 - 15. * c3;
	c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.) {
	/* 'Subroutine' F is called (kflag=0 for correct return) */
	kflag = 0;
	goto Step_F;
    }


    repeat {
	/* Step E. Exponential Sample */

	E = igraph_exp_rand(rng);/* ~ Exp(1) (standard exponential) */

	/*  sample t from the laplace 'hat'
	    (if t <= -0.6744 then pk < fk for all mu >= 10.) */
	u = 2 * igraph_rng_get_unif01(rng) - 1.;
	t = 1.8 + fsign(E, u);
	if (t > -0.6744) {
	    pois = floor(mu + s * t);
	    fk = pois;
	    difmuk = mu - fk;

	    /* 'subroutine' F is called (kflag=1 for correct return) */
	    kflag = 1;

	  Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

	    if (pois < 10) { /* use factorials from table fact[] */
		px = -mu;
		py = pow(mu, pois) / fact[(int)pois];
	    }
	    else {
		/* Case pois >= 10 uses polynomial approximation
		   a0-a7 for accuracy when advisable */
		del = one_12 / fk;
		del = del * (1. - 4.8 * del * del);
		v = difmuk / fk;
		if (fabs(v) <= 0.25)
		    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
					  v + a3) * v + a2) * v + a1) * v + a0)
			- del;
		else /* |v| > 1/4 */
		    px = fk * log(1. + v) - difmuk - del;
		py = M_1_SQRT_2PI / sqrt(fk);
	    }
	    x = (0.5 - difmuk) / s;
	    x *= x;/* x^2 */
	    fx = -0.5 * x;
	    fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
	    if (kflag > 0) {
		/* Step H. Hat acceptance (E is repeated on rejection) */
		if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
		    break;
	    } else
		/* Step Q. Quotient acceptance (rare case) */
		if (fy - u * fy <= py * exp(px - fx))
		    break;
	}/* t > -.67.. */
    }
    return pois;
}

double igraph_rgeom(igraph_rng_t *rng, double p) {
    if (ISNAN(p) || p <= 0 || p > 1) ML_ERR_return_NAN;

    return igraph_rpois(rng, igraph_exp_rand(rng) * ((1 - p) / p));
}

/* This is from nmath/rbinom.c */

#define repeat for(;;)

double igraph_rbinom(igraph_rng_t *rng, double nin, double pp)
{
    /* FIXME: These should become THREAD_specific globals : */

    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;

    static double psave = -1.0;
    static int nsave = -1;
    static int m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i,ix,k, n;

    if (!R_FINITE(nin)) ML_ERR_return_NAN;
    n = floor(nin + 0.5);
    if (n != nin) ML_ERR_return_NAN;

    if (!R_FINITE(pp) ||
	/* n=0, p=0, p=1 are not errors <TSL>*/
	n < 0 || pp < 0. || pp > 1.)	ML_ERR_return_NAN;

    if (n == 0 || pp == 0.) return 0;
    if (pp == 1.) return n;

    p = fmin(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (pp != psave || n != nsave) {
	psave = pp;
	nsave = n;
	if (np < 30.0) {
	    /* inverse cdf logic for mean less than 30 */
	    qn = pow(q, (double) n);
	    goto L_np_small;
	} else {
	    ffm = np + p;
	    m = ffm;
	    fm = m;
	    npq = np * q;
	    p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
	    xm = fm + 0.5;
	    xl = xm - p1;
	    xr = xm + p1;
	    c = 0.134 + 20.5 / (15.3 + fm);
	    al = (ffm - xl) / (ffm - xl * p);
	    xll = al * (1.0 + 0.5 * al);
	    al = (xr - ffm) / (xr * q);
	    xlr = al * (1.0 + 0.5 * al);
	    p2 = p1 * (1.0 + c + c);
	    p3 = p2 + c / xll;
	    p4 = p3 + c / xlr;
	}
    } else if (n == nsave) {
	if (np < 30.0)
	    goto L_np_small;
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
      u = igraph_rng_get_unif01(rng) * p4;
      v = igraph_rng_get_unif01(rng);
      /* triangular region */
      if (u <= p1) {
	  ix = xm - p1 * v + u;
	  goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
	  x = xl + (u - p1) / c;
	  v = v * c + 1.0 - fabs(xm - x) / p1;
	  if (v > 1.0 || v <= 0.)
	      continue;
	  ix = x;
      } else {
	  if (u > p3) {	/* right tail */
	      ix = xr - log(v) / xlr;
	      if (ix > n)
		  continue;
	      v = v * (u - p3) * xlr;
	  } else {/* left tail */
	      ix = xl + log(v) / xll;
	      if (ix < 0)
		  continue;
	      v = v * (u - p2) * xll;
	  }
      }
      /* determine appropriate way to perform accept/reject test */
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
	  /* explicit evaluation */
	  f = 1.0;
	  if (m < ix) {
	      for (i = m + 1; i <= ix; i++)
		  f *= (g / i - r);
	  } else if (m != ix) {
	      for (i = ix + 1; i <= m; i++)
		  f /= (g / i - r);
	  }
	  if (v <= f)
	      goto finis;
      } else {
	  /* squeezing using upper and lower bounds on log(f(x)) */
	  amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
	  ynorm = -k * k / (2.0 * npq);
	  alv = log(v);
	  if (alv < ynorm - amaxp)
	      goto finis;
	  if (alv <= ynorm + amaxp) {
	      /* stirling's formula to machine accuracy */
	      /* for the final acceptance/rejection test */
	      x1 = ix + 1;
	      f1 = fm + 1.0;
	      z = n + 1 - fm;
	      w = n - ix + 1.0;
	      z2 = z * z;
	      x2 = x1 * x1;
	      f2 = f1 * f1;
	      w2 = w * w;
	      if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
		  goto finis;
	  }
      }
  }

 L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

  repeat {
     ix = 0;
     f = qn;
     u = igraph_rng_get_unif01(rng);
     repeat {
	 if (u < f)
	     goto finis;
	 if (ix > 110)
	     break;
	 u -= f;
	 ix++;
	 f *= (g / ix - r);
     }
  }
 finis:
    if (psave > 0.5)
	 ix = n - ix;
  return (double)ix;
}

#endif

/**********************************************************
 * Testing purposes                                       *
 *********************************************************/

/* int main() { */

/*   int i; */

/*   RNG_BEGIN(); */

/*   for (i=0; i<1000; i++) { */
/*     printf("%li ", RNG_INTEGER(1,10)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<1000; i++) { */
/*     printf("%f ", RNG_UNIF(0,1)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<1000; i++) { */
/*     printf("%f ", RNG_NORMAL(0,5)); */
/*   } */
/*   printf("\n"); */

/*   RNG_END(); */

/*   return 0; */
/* } */
