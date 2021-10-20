/* qmd.c */

#include "env.h"
#include "qmd.h"

void genqmd(int *neqns, int xadj[], int adjncy[], int perm[],
      int invp[], int deg[], int marker[], int rchset[], int nbrhd[],
      int qsize[], int qlink[], int *nofsub)
{     static const char func[] = "genqmd";
      xassert(neqns == neqns);
      xassert(xadj == xadj);
      xassert(adjncy == adjncy);
      xassert(perm == perm);
      xassert(invp == invp);
      xassert(deg == deg);
      xassert(marker == marker);
      xassert(rchset == rchset);
      xassert(nbrhd == nbrhd);
      xassert(qsize == qsize);
      xassert(qlink == qlink);
      xassert(nofsub == nofsub);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
}

void qmdrch(int *root, int xadj[], int adjncy[], int deg[],
      int marker[], int *rchsze, int rchset[], int *nhdsze,
      int nbrhd[])
{     static const char func[] = "qmdrch";
      xassert(root == root);
      xassert(xadj == xadj);
      xassert(adjncy == adjncy);
      xassert(deg == deg);
      xassert(marker == marker);
      xassert(rchsze == rchsze);
      xassert(rchset == rchset);
      xassert(nhdsze == nhdsze);
      xassert(nbrhd == nbrhd);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
}

void qmdqt(int *root, int xadj[], int adjncy[], int marker[],
      int *rchsze, int rchset[], int nbrhd[])
{     static const char func[] = "qmdqt";
      xassert(root == root);
      xassert(xadj == xadj);
      xassert(adjncy == adjncy);
      xassert(marker == marker);
      xassert(rchsze == rchsze);
      xassert(rchset == rchset);
      xassert(nbrhd == nbrhd);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
}

void qmdupd(int xadj[], int adjncy[], int *nlist, int list[],
      int deg[], int qsize[], int qlink[], int marker[], int rchset[],
      int nbrhd[])
{     static const char func[] = "qmdupd";
      xassert(xadj == xadj);
      xassert(adjncy == adjncy);
      xassert(nlist == nlist);
      xassert(list == list);
      xassert(deg == deg);
      xassert(qsize == qsize);
      xassert(qlink == qlink);
      xassert(marker == marker);
      xassert(rchset == rchset);
      xassert(nbrhd == nbrhd);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
}

void qmdmrg(int xadj[], int adjncy[], int deg[], int qsize[],
      int qlink[], int marker[], int *deg0, int *nhdsze, int nbrhd[],
      int rchset[], int ovrlp[])
{     static const char func[] = "qmdmrg";
      xassert(xadj == xadj);
      xassert(adjncy == adjncy);
      xassert(deg == deg);
      xassert(qsize == qsize);
      xassert(qlink == qlink);
      xassert(marker == marker);
      xassert(deg0 == deg0);
      xassert(nhdsze == nhdsze);
      xassert(nbrhd == nbrhd);
      xassert(rchset == rchset);
      xassert(ovrlp == ovrlp);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
}

/* eof */
