/* glplpx03.c (OPB format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Author: Oscar Gustafsson <oscarg@isy.liu.se>.
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#define _GLPSTD_ERRNO
#define _GLPSTD_STDIO
#include "glpapi.h"
#if 0 /* 24/XII-2009; by mao */
#include "glpipp.h"
#endif

/*----------------------------------------------------------------------
-- lpx_write_pb - write problem data in (normalized) OPB format.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_write_pb(LPX *lp, const char *fname, int normalized,
--    int binarize);
--
-- *Description*
--
-- The routine lpx_write_pb writes problem data in OPB format
-- to an output text file whose name is the character string fname.
-- If normalized is non-zero the output will be generated in a
-- normalized form with sequentially numbered variables, x1, x2 etc.
-- If binarize, any integer variable will be repalzec by binary ones,
-- see ipp_binarize
--
-- *Returns*
--
-- If the operation was successful, the routine returns zero. Otherwise
-- the routine prints an error message and returns non-zero. */

#if 1 /* 24/XII-2009; by mao (disabled, because IPP was removed) */
int lpx_write_pb(LPX *lp, const char *fname, int normalized,
      int binarize)
{     xassert(lp == lp);
      xassert(fname == fname);
      xassert(normalized == normalized);
      xassert(binarize == binarize);
      xprintf("lpx_write_pb: sorry, currently this operation is not ava"
         "ilable\n");
      return 1;
}
#else
int lpx_write_pb(LPX *lp, const char *fname, int normalized,
      int binarize)
{
  FILE* fp;
  int m,n,i,j,k,o,nonfree=0, obj_dir, dbl, *ndx, row_type, emptylhs=0;
  double coeff, *val, bound, constant/*=0.0*/;
  char* objconstname = "dummy_one";
  char* emptylhsname = "dummy_zero";

  /* Variables needed for possible binarization */
  /*LPX* tlp;*/
  IPP *ipp = NULL;
  /*tlp=lp;*/

  if(binarize) /* Transform integer variables to binary ones */
    {
      ipp = ipp_create_wksp();
      ipp_load_orig(ipp, lp);
      ipp_binarize(ipp);
      lp = ipp_build_prob(ipp);
    }
  fp = fopen(fname, "w");

  if(fp!= NULL)
    {
      xprintf(
          "lpx_write_pb: writing problem in %sOPB format to `%s'...\n",
              (normalized?"normalized ":""), fname);

      m = glp_get_num_rows(lp);
      n = glp_get_num_cols(lp);
      for(i=1;i<=m;i++)
        {
          switch(glp_get_row_type(lp,i))
            {
            case GLP_LO:
            case GLP_UP:
            case GLP_FX:
              {
                nonfree += 1;
                break;
              }
            case GLP_DB:
              {
                nonfree += 2;
                break;
              }
            }
        }
      constant=glp_get_obj_coef(lp,0);
      fprintf(fp,"* #variables = %d #constraints = %d\n",
         n + (constant == 0?1:0), nonfree + (constant == 0?1:0));
      /* Objective function */
      obj_dir = glp_get_obj_dir(lp);
      fprintf(fp,"min: ");
      for(i=1;i<=n;i++)
        {
          coeff = glp_get_obj_coef(lp,i);
          if(coeff != 0.0)
            {
              if(obj_dir == GLP_MAX)
                coeff=-coeff;
              if(normalized)
                fprintf(fp, " %d x%d", (int)coeff, i);
              else
                fprintf(fp, " %d*%s", (int)coeff,
                  glp_get_col_name(lp,i));

            }
        }
      if(constant)
        {
          if(normalized)
            fprintf(fp, " %d x%d", (int)constant, n+1);
          else
            fprintf(fp, " %d*%s", (int)constant, objconstname);
        }
      fprintf(fp,";\n");

      if(normalized && !binarize)  /* Name substitution */
        {
          fprintf(fp,"* Variable name substitution:\n");
          for(j=1;j<=n;j++)
            {
              fprintf(fp, "* x%d = %s\n", j, glp_get_col_name(lp,j));
            }
          if(constant)
            fprintf(fp, "* x%d = %s\n", n+1, objconstname);
        }

      ndx = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));

      /* Constraints */
      for(j=1;j<=m;j++)
        {
          row_type=glp_get_row_type(lp,j);
          if(row_type!=GLP_FR)
            {
              if(row_type == GLP_DB)
                {
                  dbl=2;
                  row_type = GLP_UP;
                }
              else
                {
                  dbl=1;
                }
              k=glp_get_mat_row(lp, j, ndx, val);
              for(o=1;o<=dbl;o++)
                {
                  if(o==2)
                    {
                      row_type = GLP_LO;
                    }
                  if(k==0) /* Empty LHS */
                    {
                      emptylhs = 1;
                      if(normalized)
                        {
                          fprintf(fp, "0 x%d ", n+2);
                        }
                      else
                        {
                          fprintf(fp, "0*%s ", emptylhsname);
                        }
                    }

                  for(i=1;i<=k;i++)
                    {
                      if(val[i] != 0.0)
                        {

                          if(normalized)
                            {
                              fprintf(fp, "%d x%d ",
              (row_type==GLP_UP)?(-(int)val[i]):((int)val[i]), ndx[i]);
                            }
                          else
                            {
                              fprintf(fp, "%d*%s ", (int)val[i],
                                      glp_get_col_name(lp,ndx[i]));
                            }
                        }
                    }
                  switch(row_type)
                    {
                    case GLP_LO:
                      {
                        fprintf(fp, ">=");
                        bound = glp_get_row_lb(lp,j);
                        break;
                      }
                    case GLP_UP:
                      {
                        if(normalized)
                          {
                            fprintf(fp, ">=");
                            bound = -glp_get_row_ub(lp,j);
                          }
                        else
                          {
                            fprintf(fp, "<=");
                            bound = glp_get_row_ub(lp,j);
                          }

                        break;
                      }
                    case GLP_FX:
                      {
                        fprintf(fp, "=");
                        bound = glp_get_row_lb(lp,j);
                        break;
                      }
                    }
                  fprintf(fp," %d;\n",(int)bound);
                }
            }
        }
      xfree(ndx);
      xfree(val);

      if(constant)
        {
          xprintf(
        "lpx_write_pb: adding constant objective function variable\n");

          if(normalized)
            fprintf(fp, "1 x%d = 1;\n", n+1);
          else
            fprintf(fp, "1*%s = 1;\n", objconstname);
        }
      if(emptylhs)
        {
          xprintf(
            "lpx_write_pb: adding dummy variable for empty left-hand si"
            "de constraint\n");

          if(normalized)
            fprintf(fp, "1 x%d = 0;\n", n+2);
          else
            fprintf(fp, "1*%s = 0;\n", emptylhsname);
        }

    }
  else
    {
      xprintf("Problems opening file for writing: %s\n", fname);
      return(1);
    }
  fflush(fp);
  if (ferror(fp))
    {  xprintf("lpx_write_pb: can't write to `%s' - %s\n", fname,
               strerror(errno));
    goto fail;
    }
  fclose(fp);


  if(binarize)
    {
      /* delete the resultant problem object */
      if (lp != NULL) lpx_delete_prob(lp);
      /* delete MIP presolver workspace */
      if (ipp != NULL) ipp_delete_wksp(ipp);
      /*lp=tlp;*/
    }
  return 0;
 fail: if (fp != NULL) fclose(fp);
  return 1;
}
#endif

/* eof */
