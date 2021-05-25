/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#define RNG_C

#ifdef RCSID
    static const char rcsid[] = "$Id: random.cpp,v 1.15 2003/05/14 03:04:45 wilder Exp wilder $";
#endif

//________________________________________________________________________
// See the header file random.h for a description of the contents of this
// file as well as references and credits.

#include "gengraph_random.h"
#include <cmath>

using namespace std;
using namespace KW_RNG;

//________________________________________________________________________
// RNG::RNOR generates normal variates with rejection.
// nfix() generates variates after rejection in RNOR.
// Despite rejection, this method is much faster than Box-Muller.

// double RNG::nfix(slong h, ulong i)
// {
//   const double r = 3.442620f;    // The starting of the right tail
//   static double x, y;

//   for(;;) {
//     x = h * wn[i];

//     // If i == 0, handle the base strip
//     if (i==0){
//       do {
//  x = -log(rand_open01()) * 0.2904764;   // .2904764 is 1/r
//  y = -log(rand_open01());
//       } while (y + y < x * x);
//       return ((h > 0) ? r + x : -r - x);
//     }

//     // If i > 0, handle the wedges of other strips
//     if (fn[i] + rand_open01() * (fn[i - 1] - fn[i]) < exp(-.5 * x * x) )
//       return x;

//     // start all over
//     h = rand_int32();
//     i = h & 127;
//     if ((ulong) abs((sint) h) < kn[i])
//       return (h * wn[i]);
//   }

// } // RNG::nfix

// // __________________________________________________________________________
// // RNG::RNOR generates exponential variates with rejection.
// // efix() generates variates after rejection in REXP.

// double RNG::efix(ulong j, ulong i)
// {
//   double x;
//   for (;;)
//   {
//     if (i == 0)
//       return (7.69711 - log(rand_open01()));

//     x = j * we[i];
//     if (fe[i] + rand_open01() * (fe[i - 1] - fe[i]) < exp(-x))
//       return (x);

//     j = rand_int32();
//     i = (j & 255);
//     if (j < ke[i])
//       return (j * we[i]);
//   }

// } // RNG::efix

// // __________________________________________________________________________
// // This procedure creates the tables used by RNOR and REXP

// void RNG::zigset()
// {
//   const double m1 = 2147483648.0; // 2^31
//   const double m2 = 4294967296.0; // 2^32

//   const double vn = 9.91256303526217e-3;
//   const double ve = 3.949659822581572e-3;

//   double dn = 3.442619855899, tn = dn;
//   double de = 7.697117470131487, te = de;

//   int i;

//   // Set up tables for RNOR
//   double q = vn / exp(-.5 * dn * dn);
//   kn[0] = (ulong) ((dn / q) * m1);
//   kn[1] = 0;
//   wn[0] = q / m1;
//   wn[127] = dn / m1;
//   fn[0]=1.;
//   fn[127] = exp(-.5 * dn * dn);
//   for(i = 126; i >= 1; i--)
//   {
//     dn = sqrt(-2 * log(vn / dn + exp(-.5 * dn * dn)));
//     kn[i + 1] = (ulong) ((dn / tn) * m1);
//     tn = dn;
//     fn[i] = exp(-.5 * dn * dn);
//     wn[i] = dn / m1;
//   }

//   // Set up tables for REXP
//   q = ve / exp(-de);
//   ke[0] = (ulong) ((de / q) * m2);
//   ke[1] = 0;
//   we[0] = q / m2;
//   we[255] = de / m2;
//   fe[0] = 1.;
//   fe[255] = exp(-de);
//   for (i = 254; i >= 1; i--)
//   {
//     de = -log(ve / de + exp(-de));
//     ke[i+1] = (ulong) ((de / te) * m2);
//     te = de;
//     fe[i] = exp(-de);
//     we[i] = de / m2;
//   }

// } // RNG::zigset

// // __________________________________________________________________________
// // Generate a gamma variate with parameters 'shape' and 'scale'

// double RNG::gamma(double shape, double scale)
// {
//   if (shape < 1)
//     return gamma(shape + 1, scale) * pow(rand_open01(), 1.0 / shape);

//   const double d = shape - 1.0 / 3.0;
//   const double c = 1.0 / sqrt(9.0 * d);
//   double x, v, u;
//   for (;;) {
//     do {
//       x = RNOR();
//       v = 1.0 + c * x;
//     } while (v <= 0.0);
//     v = v * v * v;
//     u = rand_open01();
//     if (u < 1.0 - 0.0331 * x * x * x * x)
//       return (d * v / scale);
//     if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v)))
//       return (d * v / scale);
//   }

// } // RNG::gamma

// // __________________________________________________________________________
// // gammalog returns the logarithm of the gamma function.  From Numerical
// // Recipes.

// double gammalog(double xx)
// {
//   static double cof[6]={
//     76.18009172947146, -86.50532032941677, 24.01409824083091,
//     -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

//   double x = xx;
//   double y = xx;
//   double tmp = x + 5.5;
//   tmp -= (x + 0.5) * log(tmp);
//   double ser=1.000000000190015;
//   for (int j=0; j<=5; j++)
//     ser += cof[j] / ++y;
//   return -tmp + log(2.5066282746310005 * ser / x);
// }

// // __________________________________________________________________________
// // Generate a Poisson variate
// // This is essentially the algorithm from Numerical Recipes

// double RNG::poisson(double lambda)
// {
//   static double sq, alxm, g, oldm = -1.0;
//   double em, t, y;

//   if (lambda < 12.0) {
//     if (lambda != oldm) {
//       oldm = lambda;
//       g = exp(-lambda);
//     }
//     em = -1;
//     t = 1.0;
//     do {
//       ++em;
//       t *= rand_open01();
//     } while (t > g);
//   } else {
//     if (lambda != oldm) {
//       oldm = lambda;
//       sq = sqrt(2.0 * lambda);
//       alxm = log(lambda);
//       g = lambda * alxm - gammalog(lambda + 1.0);
//     }
//     do {
//       do {
//  y = tan(PI * rand_open01());
//  em = sq * y + lambda;
//       } while (em < 0.0);
//       em = floor(em);
//       t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammalog(em + 1.0)-g);
//     } while (rand_open01() > t);
//   }
//   return em;

// } // RNG::poisson

// // __________________________________________________________________________
// // Generate a binomial variate
// // This is essentially the algorithm from Numerical Recipes

// int RNG::binomial(double pp, int n)
// {
//   if(n==0) return 0;
//   if(pp==0.0) return 0;
//   if(pp==1.0) return n;
//   double p = (pp<0.5 ? pp : 1.0-pp);
//   double am = n*p;
//   int bnl = 0;
//   if(n<25) {
//     for(int j=n; j--; ) if(rand_closed01()<p) ++bnl;
//   }
//   else if(am<1.0) {
//     double g = exp(-am);
//     double t = 1.0;
//     for (; bnl<n; bnl++) if((t*=rand_closed01())<g) break;
//   }
//   else {
//     double en = n;
//     double oldg = gammalog(en + 1.0);
//     double pc = 1.0 - p;
//     double sq = sqrt(2.0 * am * pc);
//     double y, em, t;
//     do {
//       do {
//         double angle = PI * rand_halfclosed01();
//          y = tan(angle);
//         em = sq * y + am;
//       } while (em < 0.0 || em >= en + 1.0);
//       em = floor(em);
//       t = 1.2 * sq * (1 + y * y) * exp(oldg - gammalog(em + 1.0) -
//           gammalog(en - em + 1.0) + em * log(p) + (en - em) * log(pc));
//     } while (rand_closed01() > t);
//     bnl = int(em);
//   }
//   if (p!=pp) bnl=n-bnl;
//   return bnl;
// } // RNG::binomial

// __________________________________________________________________________
// rng.C
