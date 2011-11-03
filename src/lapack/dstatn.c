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

#include "f2c.h"


/*     %---------------------------------------------%   
       | Initialize statistic and timing information |   
       | for nonsymmetric Arnoldi code.              |   
       %---------------------------------------------%   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \SCCS Information: @(#)   
   FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2   

   Subroutine */ int igraphdstatn_(void)
{
    integer nbx, nopx;
    real trvec, tmvbx, tnaup2, tgetv0, tneigh;
    integer nitref;
    real tnaupd, titref, tnaitr, tngets, tnapps, tnconv;
    integer nrorth, nrstrt;
    real tmvopx;


/*     %--------------------------------%   
       | See stat.doc for documentation |   
       %--------------------------------%   


       %-----------------------%   
       | Executable Statements |   
       %-----------------------% */

    nopx = 0;
    nbx = 0;
    nrorth = 0;
    nitref = 0;
    nrstrt = 0;

    tnaupd = 0.f;
    tnaup2 = 0.f;
    tnaitr = 0.f;
    tneigh = 0.f;
    tngets = 0.f;
    tnapps = 0.f;
    tnconv = 0.f;
    titref = 0.f;
    tgetv0 = 0.f;
    trvec = 0.f;

/*     %----------------------------------------------------%   
       | User time including reverse communication overhead |   
       %----------------------------------------------------% */

    tmvopx = 0.f;
    tmvbx = 0.f;

    return 0;


/*     %---------------%   
       | End of dstatn |   
       %---------------% */

} /* igraphdstatn_ */

