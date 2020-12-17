/* glpmpl06.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#pragma clang diagnostic ignored "-Wself-assign"
#endif

#define _GLPSTD_ERRNO
#define _GLPSTD_STDIO
#include "glpmpl.h"
#include "glpsql.h"

/**********************************************************************/

#define CSV_FIELD_MAX 50
/* maximal number of fields in record */

#define CSV_FDLEN_MAX 100
/* maximal field length */

struct csv
{     /* comma-separated values file */
      int mode;
      /* 'R' = reading; 'W' = writing */
      char *fname;
      /* name of csv file */
      FILE *fp;
      /* stream assigned to csv file */
      jmp_buf jump;
      /* address for non-local go to in case of error */
      int count;
      /* record count */
      /*--------------------------------------------------------------*/
      /* used only for input csv file */
      int c;
      /* current character or EOF */
      int what;
      /* current marker: */
#define CSV_EOF   0  /* end-of-file */
#define CSV_EOR   1  /* end-of-record */
#define CSV_NUM   2  /* floating-point number */
#define CSV_STR   3  /* character string */
      char field[CSV_FDLEN_MAX+1];
      /* current field just read */
      int nf;
      /* number of fields in the csv file */
      int ref[1+CSV_FIELD_MAX];
      /* ref[k] = k', if k-th field of the csv file corresponds to
         k'-th field in the table statement; if ref[k] = 0, k-th field
         of the csv file is ignored */
#if 1 /* 01/VI-2010 */
      int nskip;
      /* number of comment records preceding the header record */
#endif
};

#undef read_char

static void read_char(struct csv *csv)
{     /* read character from csv data file */
      int c;
      xassert(csv->c != EOF);
      if (csv->c == '\n') csv->count++;
loop: c = fgetc(csv->fp);
      if (ferror(csv->fp))
      {  xprintf("%s:%d: read error - %s\n", csv->fname, csv->count,
            strerror(errno));
         longjmp(csv->jump, 0);
      }
      if (feof(csv->fp))
      {  if (csv->c == '\n')
         {  csv->count--;
            c = EOF;
         }
         else
         {  xprintf("%s:%d: warning: missing final end-of-line\n",
               csv->fname, csv->count);
            c = '\n';
         }
      }
      else if (c == '\r')
         goto loop;
      else if (c == '\n')
         ;
      else if (iscntrl(c))
      {  xprintf("%s:%d: invalid control character 0x%02X\n",
            csv->fname, csv->count, c);
         longjmp(csv->jump, 0);
      }
      csv->c = c;
      return;
}

static void read_field(struct csv *csv)
{     /* read field from csv data file */
      /* check for end of file */
      if (csv->c == EOF)
      {  csv->what = CSV_EOF;
         strcpy(csv->field, "EOF");
         goto done;
      }
      /* check for end of record */
      if (csv->c == '\n')
      {  csv->what = CSV_EOR;
         strcpy(csv->field, "EOR");
         read_char(csv);
         if (csv->c == ',')
err1:    {  xprintf("%s:%d: empty field not allowed\n", csv->fname,
               csv->count);
            longjmp(csv->jump, 0);
         }
         if (csv->c == '\n')
         {  xprintf("%s:%d: empty record not allowed\n", csv->fname,
               csv->count);
            longjmp(csv->jump, 0);
         }
#if 1 /* 01/VI-2010 */
         /* skip comment records; may appear only before the very first
            record containing field names */
         if (csv->c == '#' && csv->count == 1)
         {  while (csv->c == '#')
            {  while (csv->c != '\n')
                  read_char(csv);
               read_char(csv);
               csv->nskip++;
            }
         }
#endif
         goto done;
      }
      /* skip comma before next field */
      if (csv->c == ',')
         read_char(csv);
      /* read field */
      if (csv->c == '\'' || csv->c == '"')
      {  /* read a field enclosed in quotes */
         int quote = csv->c, len = 0;
         csv->what = CSV_STR;
         /* skip opening quote */
         read_char(csv);
         /* read field characters within quotes */
         for (;;)
         {  /* check for closing quote and read it */
            if (csv->c == quote)
            {  read_char(csv);
               if (csv->c == quote)
                  ;
               else if (csv->c == ',' || csv->c == '\n')
                  break;
               else
               {  xprintf("%s:%d: invalid field\n", csv->fname,
                     csv->count);
                  longjmp(csv->jump, 0);
               }
            }
            /* check the current field length */
            if (len == CSV_FDLEN_MAX)
err2:       {  xprintf("%s:%d: field too long\n", csv->fname,
                  csv->count);
               longjmp(csv->jump, 0);
            }
            /* add the current character to the field */
            csv->field[len++] = (char)csv->c;
            /* read the next character */
            read_char(csv);
         }
         /* the field has been read */
         if (len == 0) goto err1;
         csv->field[len] = '\0';
      }
      else
      {  /* read a field not enclosed in quotes */
         int len = 0;
         double temp;
         csv->what = CSV_NUM;
         while (!(csv->c == ',' || csv->c == '\n'))
         {  /* quotes within the field are not allowed */
            if (csv->c == '\'' || csv->c == '"')
            {  xprintf("%s:%d: invalid use of single or double quote wi"
                  "thin field\n", csv->fname, csv->count);
               longjmp(csv->jump, 0);
            }
            /* check the current field length */
            if (len == CSV_FDLEN_MAX) goto err2;
            /* add the current character to the field */
            csv->field[len++] = (char)csv->c;
            /* read the next character */
            read_char(csv);
         }
         /* the field has been read */
         if (len == 0) goto err1;
         csv->field[len] = '\0';
         /* check the field type */
         if (str2num(csv->field, &temp)) csv->what = CSV_STR;
      }
done: return;
}

static struct csv *csv_open_file(TABDCA *dca, int mode)
{     /* open csv data file */
      struct csv *csv;
      /* create control structure */
      csv = xmalloc(sizeof(struct csv));
      csv->mode = mode;
      csv->fname = NULL;
      csv->fp = NULL;
      if (setjmp(csv->jump)) goto fail;
      csv->count = 0;
      csv->c = '\n';
      csv->what = 0;
      csv->field[0] = '\0';
      csv->nf = 0;
      /* try to open the csv data file */
      if (mpl_tab_num_args(dca) < 2)
      {  xprintf("csv_driver: file name not specified\n");
         longjmp(csv->jump, 0);
      }
      csv->fname = xmalloc(strlen(mpl_tab_get_arg(dca, 2))+1);
      strcpy(csv->fname, mpl_tab_get_arg(dca, 2));
      if (mode == 'R')
      {  /* open the file for reading */
         int k;
         csv->fp = fopen(csv->fname, "r");
         if (csv->fp == NULL)
         {  xprintf("csv_driver: unable to open %s - %s\n",
               csv->fname, strerror(errno));
            longjmp(csv->jump, 0);
         }
#if 1 /* 01/VI-2010 */
         csv->nskip = 0;
#endif
         /* skip fake new-line */
         read_field(csv);
         xassert(csv->what == CSV_EOR);
         /* read field names */
         xassert(csv->nf == 0);
         for (;;)
         {  read_field(csv);
            if (csv->what == CSV_EOR)
               break;
            if (csv->what != CSV_STR)
            {  xprintf("%s:%d: invalid field name\n", csv->fname,
                  csv->count);
               longjmp(csv->jump, 0);
            }
            if (csv->nf == CSV_FIELD_MAX)
            {  xprintf("%s:%d: too many fields\n", csv->fname,
                  csv->count);
               longjmp(csv->jump, 0);
            }
            csv->nf++;
            /* find corresponding field in the table statement */
            for (k = mpl_tab_num_flds(dca); k >= 1; k--)
            {  if (strcmp(mpl_tab_get_name(dca, k), csv->field) == 0)
                  break;
            }
            csv->ref[csv->nf] = k;
         }
         /* find dummy RECNO field in the table statement */
         for (k = mpl_tab_num_flds(dca); k >= 1; k--)
            if (strcmp(mpl_tab_get_name(dca, k), "RECNO") == 0) break;
         csv->ref[0] = k;
      }
      else if (mode == 'W')
      {  /* open the file for writing */
         int k, nf;
         csv->fp = fopen(csv->fname, "w");
         if (csv->fp == NULL)
         {  xprintf("csv_driver: unable to create %s - %s\n",
               csv->fname, strerror(errno));
            longjmp(csv->jump, 0);
         }
         /* write field names */
         nf = mpl_tab_num_flds(dca);
         for (k = 1; k <= nf; k++)
            fprintf(csv->fp, "%s%c", mpl_tab_get_name(dca, k),
               k < nf ? ',' : '\n');
         csv->count++;
      }
      else
         xassert(mode != mode);
      /* the file has been open */
      return csv;
fail: /* the file cannot be open */
      if (csv->fname != NULL) xfree(csv->fname);
      if (csv->fp != NULL) fclose(csv->fp);
      xfree(csv);
      return NULL;
}

static int csv_read_record(TABDCA *dca, struct csv *csv)
{     /* read next record from csv data file */
      int k, ret = 0;
      xassert(csv->mode == 'R');
      if (setjmp(csv->jump))
      {  ret = 1;
         goto done;
      }
      /* read dummy RECNO field */
      if (csv->ref[0] > 0)
#if 0 /* 01/VI-2010 */
         mpl_tab_set_num(dca, csv->ref[0], csv->count-1);
#else
         mpl_tab_set_num(dca, csv->ref[0], csv->count-csv->nskip-1);
#endif
      /* read fields */
      for (k = 1; k <= csv->nf; k++)
      {  read_field(csv);
         if (csv->what == CSV_EOF)
         {  /* end-of-file reached */
            xassert(k == 1);
            ret = -1;
            goto done;
         }
         else if (csv->what == CSV_EOR)
         {  /* end-of-record reached */
            int lack = csv->nf - k + 1;
            if (lack == 1)
               xprintf("%s:%d: one field missing\n", csv->fname,
                  csv->count);
            else
               xprintf("%s:%d: %d fields missing\n", csv->fname,
                  csv->count, lack);
            longjmp(csv->jump, 0);
         }
         else if (csv->what == CSV_NUM)
         {  /* floating-point number */
            if (csv->ref[k] > 0)
            {  double num;
               xassert(str2num(csv->field, &num) == 0);
               mpl_tab_set_num(dca, csv->ref[k], num);
            }
         }
         else if (csv->what == CSV_STR)
         {  /* character string */
            if (csv->ref[k] > 0)
               mpl_tab_set_str(dca, csv->ref[k], csv->field);
         }
         else
            xassert(csv != csv);
      }
      /* now there must be NL */
      read_field(csv);
      xassert(csv->what != CSV_EOF);
      if (csv->what != CSV_EOR)
      {  xprintf("%s:%d: too many fields\n", csv->fname, csv->count);
         longjmp(csv->jump, 0);
      }
done: return ret;
}

static int csv_write_record(TABDCA *dca, struct csv *csv)
{     /* write next record to csv data file */
      int k, nf, ret = 0;
      const char *c;
      xassert(csv->mode == 'W');
      nf = mpl_tab_num_flds(dca);
      for (k = 1; k <= nf; k++)
      {  switch (mpl_tab_get_type(dca, k))
         {  case 'N':
               fprintf(csv->fp, "%.*g", DBL_DIG,
                  mpl_tab_get_num(dca, k));
               break;
            case 'S':
               fputc('"', csv->fp);
               for (c = mpl_tab_get_str(dca, k); *c != '\0'; c++)
               {  if (*c == '"')
                     fputc('"', csv->fp), fputc('"', csv->fp);
                  else
                     fputc(*c, csv->fp);
               }
               fputc('"', csv->fp);
               break;
            default:
               xassert(dca != dca);
         }
         fputc(k < nf ? ',' : '\n', csv->fp);
      }
      csv->count++;
      if (ferror(csv->fp))
      {  xprintf("%s:%d: write error - %s\n", csv->fname, csv->count,
            strerror(errno));
         ret = 1;
      }
      return ret;
}

static int csv_close_file(TABDCA *dca, struct csv *csv)
{     /* close csv data file */
      int ret = 0;
      xassert(dca == dca);
      if (csv->mode == 'W')
      {  fflush(csv->fp);
         if (ferror(csv->fp))
         {  xprintf("%s:%d: write error - %s\n", csv->fname,
               csv->count, strerror(errno));
            ret = 1;
         }
      }
      xfree(csv->fname);
      fclose(csv->fp);
      xfree(csv);
      return ret;
}

/**********************************************************************/

#define DBF_FIELD_MAX 50
/* maximal number of fields in record */

#define DBF_FDLEN_MAX 100
/* maximal field length */

struct dbf
{     /* xBASE data file */
      int mode;
      /* 'R' = reading; 'W' = writing */
      char *fname;
      /* name of xBASE file */
      FILE *fp;
      /* stream assigned to xBASE file */
      jmp_buf jump;
      /* address for non-local go to in case of error */
      int offset;
      /* offset of a byte to be read next */
      int count;
      /* record count */
      int nf;
      /* number of fields */
      int ref[1+DBF_FIELD_MAX];
      /* ref[k] = k', if k-th field of the csv file corresponds to
         k'-th field in the table statement; if ref[k] = 0, k-th field
         of the csv file is ignored */
      int type[1+DBF_FIELD_MAX];
      /* type[k] is type of k-th field */
      int len[1+DBF_FIELD_MAX];
      /* len[k] is length of k-th field */
      int prec[1+DBF_FIELD_MAX];
      /* prec[k] is precision of k-th field */
};

static int read_byte(struct dbf *dbf)
{     /* read byte from xBASE data file */
      int b;
      b = fgetc(dbf->fp);
      if (ferror(dbf->fp))
      {  xprintf("%s:0x%X: read error - %s\n", dbf->fname,
            dbf->offset, strerror(errno));
         longjmp(dbf->jump, 0);
      }
      if (feof(dbf->fp))
      {  xprintf("%s:0x%X: unexpected end of file\n", dbf->fname,
            dbf->offset);
         longjmp(dbf->jump, 0);
      }
      xassert(0x00 <= b && b <= 0xFF);
      dbf->offset++;
      return b;
}

static void read_header(TABDCA *dca, struct dbf *dbf)
{     /* read xBASE data file header */
      int b, j, k, recl;
      char name[10+1];
      /* (ignored) */
      for (j = 1; j <= 10; j++)
         read_byte(dbf);
      /* length of each record, in bytes */
      recl = read_byte(dbf);
      recl += read_byte(dbf) << 8;
      /* (ignored) */
      for (j = 1; j <= 20; j++)
         read_byte(dbf);
      /* field descriptor array */
      xassert(dbf->nf == 0);
      for (;;)
      {  /* check for end of array */
         b = read_byte(dbf);
         if (b == 0x0D) break;
         if (dbf->nf == DBF_FIELD_MAX)
         {  xprintf("%s:0x%X: too many fields\n", dbf->fname,
               dbf->offset);
            longjmp(dbf->jump, 0);
         }
         dbf->nf++;
         /* field name */
         name[0] = (char)b;
         for (j = 1; j < 10; j++)
         {  b = read_byte(dbf);
            name[j] = (char)b;
         }
         name[10] = '\0';
         b = read_byte(dbf);
         if (b != 0x00)
         {  xprintf("%s:0x%X: invalid field name\n", dbf->fname,
               dbf->offset);
            longjmp(dbf->jump, 0);
         }
         /* find corresponding field in the table statement */
         for (k = mpl_tab_num_flds(dca); k >= 1; k--)
            if (strcmp(mpl_tab_get_name(dca, k), name) == 0) break;
         dbf->ref[dbf->nf] = k;
         /* field type */
         b = read_byte(dbf);
         if (!(b == 'C' || b == 'N'))
         {  xprintf("%s:0x%X: invalid field type\n", dbf->fname,
               dbf->offset);
            longjmp(dbf->jump, 0);
         }
         dbf->type[dbf->nf] = b;
         /* (ignored) */
         for (j = 1; j <= 4; j++)
            read_byte(dbf);
         /* field length */
         b = read_byte(dbf);
         if (b == 0)
         {  xprintf("%s:0x%X: invalid field length\n", dbf->fname,
               dbf->offset);
            longjmp(dbf->jump, 0);
         }
         if (b > DBF_FDLEN_MAX)
         {  xprintf("%s:0x%X: field too long\n", dbf->fname,
               dbf->offset);
            longjmp(dbf->jump, 0);
         }
         dbf->len[dbf->nf] = b;
         recl -= b;
         /* (ignored) */
         for (j = 1; j <= 15; j++)
            read_byte(dbf);
      }
      if (recl != 1)
      {  xprintf("%s:0x%X: invalid file header\n", dbf->fname,
            dbf->offset);
         longjmp(dbf->jump, 0);
      }
      /* find dummy RECNO field in the table statement */
      for (k = mpl_tab_num_flds(dca); k >= 1; k--)
         if (strcmp(mpl_tab_get_name(dca, k), "RECNO") == 0) break;
      dbf->ref[0] = k;
      return;
}

static void parse_third_arg(TABDCA *dca, struct dbf *dbf)
{     /* parse xBASE file format (third argument) */
      int j, k, temp;
      const char *arg;
      dbf->nf = mpl_tab_num_flds(dca);
      arg = mpl_tab_get_arg(dca, 3), j = 0;
      for (k = 1; k <= dbf->nf; k++)
      {  /* parse specification of k-th field */
         if (arg[j] == '\0')
         {  xprintf("xBASE driver: field %s: specification missing\n",
               mpl_tab_get_name(dca, k));
            longjmp(dbf->jump, 0);
         }
         /* parse field type */
         if (arg[j] == 'C' || arg[j] == 'N')
            dbf->type[k] = arg[j], j++;
         else
         {  xprintf("xBASE driver: field %s: invalid field type\n",
               mpl_tab_get_name(dca, k));
            longjmp(dbf->jump, 0);
         }
         /* check for left parenthesis */
         if (arg[j] == '(')
            j++;
         else
err:     {  xprintf("xBASE driver: field %s: invalid field format\n",
               mpl_tab_get_name(dca, k));
            longjmp(dbf->jump, 0);
         }
         /* parse field length */
         temp = 0;
         while (isdigit(arg[j]))
         {  if (temp > DBF_FDLEN_MAX) break;
            temp = 10 * temp + (arg[j] - '0'), j++;
         }
         if (!(1 <= temp && temp <= DBF_FDLEN_MAX))
         {  xprintf("xBASE driver: field %s: invalid field length\n",
               mpl_tab_get_name(dca, k));
            longjmp(dbf->jump, 0);
         }
         dbf->len[k] = temp;
         /* parse optional field precision */
         if (dbf->type[k] == 'N' && arg[j] == ',')
         {  j++;
            temp = 0;
            while (isdigit(arg[j]))
            {  if (temp > dbf->len[k]) break;
               temp = 10 * temp + (arg[j] - '0'), j++;
            }
            if (temp > dbf->len[k])
            {  xprintf("xBASE driver: field %s: invalid field precision"
                  "\n", mpl_tab_get_name(dca, k));
               longjmp(dbf->jump, 0);
            }
            dbf->prec[k] = temp;
         }
         else
            dbf->prec[k] = 0;
         /* check for right parenthesis */
         if (arg[j] == ')')
            j++;
         else
            goto err;
      }
      /* ignore other specifications */
      return;
}

static void write_byte(struct dbf *dbf, int b)
{     /* write byte to xBASE data file */
      fputc(b, dbf->fp);
      dbf->offset++;
      return;
}

static void write_header(TABDCA *dca, struct dbf *dbf)
{     /* write xBASE data file header */
      int j, k, temp;
      const char *name;
      /* version number */
      write_byte(dbf, 0x03 /* file without DBT */);
      /* date of last update (YYMMDD) */
      write_byte(dbf, 70 /* 1970 */);
      write_byte(dbf, 1 /* January */);
      write_byte(dbf, 1 /* 1st */);
      /* number of records (unknown so far) */
      for (j = 1; j <= 4; j++)
         write_byte(dbf, 0xFF);
      /* length of the header, in bytes */
      temp = 32 + dbf->nf * 32 + 1;
      write_byte(dbf, temp);
      write_byte(dbf, temp >> 8);
      /* length of each record, in bytes */
      temp = 1;
      for (k = 1; k <= dbf->nf; k++)
         temp += dbf->len[k];
      write_byte(dbf, temp);
      write_byte(dbf, temp >> 8);
      /* (reserved) */
      for (j = 1; j <= 20; j++)
         write_byte(dbf, 0x00);
      /* field descriptor array */
      for (k = 1; k <= dbf->nf; k++)
      {  /* field name (terminated by 0x00) */
         name = mpl_tab_get_name(dca, k);
         for (j = 0; j < 10 && name[j] != '\0'; j++)
            write_byte(dbf, name[j]);
         for (j = j; j < 11; j++)
            write_byte(dbf, 0x00);
         /* field type */
         write_byte(dbf, dbf->type[k]);
         /* (reserved) */
         for (j = 1; j <= 4; j++)
            write_byte(dbf, 0x00);
         /* field length */
         write_byte(dbf, dbf->len[k]);
         /* field precision */
         write_byte(dbf, dbf->prec[k]);
         /* (reserved) */
         for (j = 1; j <= 14; j++)
            write_byte(dbf, 0x00);
      }
      /* end of header */
      write_byte(dbf, 0x0D);
      return;
}

static struct dbf *dbf_open_file(TABDCA *dca, int mode)
{     /* open xBASE data file */
      struct dbf *dbf;
      /* create control structure */
      dbf = xmalloc(sizeof(struct dbf));
      dbf->mode = mode;
      dbf->fname = NULL;
      dbf->fp = NULL;
      if (setjmp(dbf->jump)) goto fail;
      dbf->offset = 0;
      dbf->count = 0;
      dbf->nf = 0;
      /* try to open the xBASE data file */
      if (mpl_tab_num_args(dca) < 2)
      {  xprintf("xBASE driver: file name not specified\n");
         longjmp(dbf->jump, 0);
      }
      dbf->fname = xmalloc(strlen(mpl_tab_get_arg(dca, 2))+1);
      strcpy(dbf->fname, mpl_tab_get_arg(dca, 2));
      if (mode == 'R')
      {  /* open the file for reading */
         dbf->fp = fopen(dbf->fname, "rb");
         if (dbf->fp == NULL)
         {  xprintf("xBASE driver: unable to open %s - %s\n",
               dbf->fname, strerror(errno));
            longjmp(dbf->jump, 0);
         }
         read_header(dca, dbf);
      }
      else if (mode == 'W')
      {  /* open the file for writing */
         if (mpl_tab_num_args(dca) < 3)
         {  xprintf("xBASE driver: file format not specified\n");
            longjmp(dbf->jump, 0);
         }
         parse_third_arg(dca, dbf);
         dbf->fp = fopen(dbf->fname, "wb");
         if (dbf->fp == NULL)
         {  xprintf("xBASE driver: unable to create %s - %s\n",
               dbf->fname, strerror(errno));
            longjmp(dbf->jump, 0);
         }
         write_header(dca, dbf);
      }
      else
         xassert(mode != mode);
      /* the file has been open */
      return dbf;
fail: /* the file cannot be open */
      if (dbf->fname != NULL) xfree(dbf->fname);
      if (dbf->fp != NULL) fclose(dbf->fp);
      xfree(dbf);
      return NULL;
}

static int dbf_read_record(TABDCA *dca, struct dbf *dbf)
{     /* read next record from xBASE data file */
      int b, j, k, ret = 0;
      char buf[DBF_FDLEN_MAX+1];
      xassert(dbf->mode == 'R');
      if (setjmp(dbf->jump))
      {  ret = 1;
         goto done;
      }
      /* check record flag */
      b = read_byte(dbf);
      if (b == 0x1A)
      {  /* end of data */
         ret = -1;
         goto done;
      }
      if (b != 0x20)
      {  xprintf("%s:0x%X: invalid record flag\n", dbf->fname,
            dbf->offset);
         longjmp(dbf->jump, 0);
      }
      /* read dummy RECNO field */
      if (dbf->ref[0] > 0)
         mpl_tab_set_num(dca, dbf->ref[0], dbf->count+1);
      /* read fields */
      for (k = 1; k <= dbf->nf; k++)
      {  /* read k-th field */
         for (j = 0; j < dbf->len[k]; j++)
            buf[j] = (char)read_byte(dbf);
         buf[dbf->len[k]] = '\0';
         /* set field value */
         if (dbf->type[k] == 'C')
         {  /* character field */
            if (dbf->ref[k] > 0)
               mpl_tab_set_str(dca, dbf->ref[k], strtrim(buf));
         }
         else if (dbf->type[k] == 'N')
         {  /* numeric field */
            if (dbf->ref[k] > 0)
            {  double num;
               strspx(buf);
               xassert(str2num(buf, &num) == 0);
               mpl_tab_set_num(dca, dbf->ref[k], num);
            }
         }
         else
            xassert(dbf != dbf);
      }
      /* increase record count */
      dbf->count++;
done: return ret;
}

static int dbf_write_record(TABDCA *dca, struct dbf *dbf)
{     /* write next record to xBASE data file */
      int j, k, ret = 0;
      char buf[255+1];
      xassert(dbf->mode == 'W');
      if (setjmp(dbf->jump))
      {  ret = 1;
         goto done;
      }
      /* record flag */
      write_byte(dbf, 0x20);
      xassert(dbf->nf == mpl_tab_num_flds(dca));
      for (k = 1; k <= dbf->nf; k++)
      {  if (dbf->type[k] == 'C')
         {  /* character field */
            const char *str;
            if (mpl_tab_get_type(dca, k) == 'N')
            {  sprintf(buf, "%.*g", DBL_DIG, mpl_tab_get_num(dca, k));
               str = buf;
            }
            else if (mpl_tab_get_type(dca, k) == 'S')
               str = mpl_tab_get_str(dca, k);
            else
               xassert(dca != dca);
            if ((int)strlen(str) > dbf->len[k])
            {  xprintf("xBASE driver: field %s: cannot convert %.15s..."
                  " to field format\n", mpl_tab_get_name(dca, k), str);
               longjmp(dbf->jump, 0);
            }
            for (j = 0; j < dbf->len[k] && str[j] != '\0'; j++)
                write_byte(dbf, str[j]);
            for (j = j; j < dbf->len[k]; j++)
                write_byte(dbf, ' ');
         }
         else if (dbf->type[k] == 'N')
         {  /* numeric field */
            double num = mpl_tab_get_num(dca, k);
            if (fabs(num) > 1e20)
err:        {  xprintf("xBASE driver: field %s: cannot convert %g to fi"
                  "eld format\n", mpl_tab_get_name(dca, k), num);
               longjmp(dbf->jump, 0);
            }
            sprintf(buf, "%*.*f", dbf->len[k], dbf->prec[k], num);
            xassert(strlen(buf) < sizeof(buf));
            if ((int)strlen(buf) != dbf->len[k]) goto err;
            for (j = 0; j < dbf->len[k]; j++)
               write_byte(dbf, buf[j]);
         }
         else
            xassert(dbf != dbf);
      }
      /* increase record count */
      dbf->count++;
done: return ret;
}

static int dbf_close_file(TABDCA *dca, struct dbf *dbf)
{     /* close xBASE data file */
      int ret = 0;
      xassert(dca == dca);
      if (dbf->mode == 'W')
      {  if (setjmp(dbf->jump))
         {  ret = 1;
            goto skip;
         }
         /* end-of-file flag */
         write_byte(dbf, 0x1A);
         /* number of records */
         dbf->offset = 4;
         if (fseek(dbf->fp, dbf->offset, SEEK_SET))
         {  xprintf("%s:0x%X: seek error - %s\n", dbf->fname,
               dbf->offset, strerror(errno));
            longjmp(dbf->jump, 0);
         }
         write_byte(dbf, dbf->count);
         write_byte(dbf, dbf->count >> 8);
         write_byte(dbf, dbf->count >> 16);
         write_byte(dbf, dbf->count >> 24);
         fflush(dbf->fp);
         if (ferror(dbf->fp))
         {  xprintf("%s:0x%X: write error - %s\n", dbf->fname,
               dbf->offset, strerror(errno));
            longjmp(dbf->jump, 0);
         }
skip:    ;
      }
      xfree(dbf->fname);
      fclose(dbf->fp);
      xfree(dbf);
      return ret;
}

/**********************************************************************/

#define TAB_CSV   1
#define TAB_XBASE 2
#define TAB_ODBC  3
#define TAB_MYSQL 4

void mpl_tab_drv_open(MPL *mpl, int mode)
{     TABDCA *dca = mpl->dca;
      xassert(dca->id == 0);
      xassert(dca->link == NULL);
      xassert(dca->na >= 1);
      if (strcmp(dca->arg[1], "CSV") == 0)
      {  dca->id = TAB_CSV;
         dca->link = csv_open_file(dca, mode);
      }
      else if (strcmp(dca->arg[1], "xBASE") == 0)
      {  dca->id = TAB_XBASE;
         dca->link = dbf_open_file(dca, mode);
      }
      else if (strcmp(dca->arg[1], "ODBC") == 0 ||
               strcmp(dca->arg[1], "iODBC") == 0)
      {  dca->id = TAB_ODBC;
         dca->link = db_iodbc_open(dca, mode);
      }
      else if (strcmp(dca->arg[1], "MySQL") == 0)
      {  dca->id = TAB_MYSQL;
         dca->link = db_mysql_open(dca, mode);
      }
      else
         xprintf("Invalid table driver `%s'\n", dca->arg[1]);
      if (dca->link == NULL)
         error(mpl, "error on opening table %s",
            mpl->stmt->u.tab->name);
      return;
}

int mpl_tab_drv_read(MPL *mpl)
{     TABDCA *dca = mpl->dca;
      int ret;
      switch (dca->id)
      {  case TAB_CSV:
            ret = csv_read_record(dca, dca->link);
            break;
         case TAB_XBASE:
            ret = dbf_read_record(dca, dca->link);
            break;
         case TAB_ODBC:
            ret = db_iodbc_read(dca, dca->link);
            break;
         case TAB_MYSQL:
            ret = db_mysql_read(dca, dca->link);
            break;
         default:
            xassert(dca != dca);
      }
      if (ret > 0)
         error(mpl, "error on reading data from table %s",
            mpl->stmt->u.tab->name);
      return ret;
}

void mpl_tab_drv_write(MPL *mpl)
{     TABDCA *dca = mpl->dca;
      int ret;
      switch (dca->id)
      {  case TAB_CSV:
            ret = csv_write_record(dca, dca->link);
            break;
         case TAB_XBASE:
            ret = dbf_write_record(dca, dca->link);
            break;
         case TAB_ODBC:
            ret = db_iodbc_write(dca, dca->link);
            break;
         case TAB_MYSQL:
            ret = db_mysql_write(dca, dca->link);
            break;
         default:
            xassert(dca != dca);
      }
      if (ret)
         error(mpl, "error on writing data to table %s",
            mpl->stmt->u.tab->name);
      return;
}

void mpl_tab_drv_close(MPL *mpl)
{     TABDCA *dca = mpl->dca;
      int ret;
      switch (dca->id)
      {  case TAB_CSV:
            ret = csv_close_file(dca, dca->link);
            break;
         case TAB_XBASE:
            ret = dbf_close_file(dca, dca->link);
            break;
         case TAB_ODBC:
            ret = db_iodbc_close(dca, dca->link);
            break;
         case TAB_MYSQL:
            ret = db_mysql_close(dca, dca->link);
            break;
         default:
            xassert(dca != dca);
      }
      dca->id = 0;
      dca->link = NULL;
      if (ret)
         error(mpl, "error on closing table %s",
            mpl->stmt->u.tab->name);
      return;
}

/* eof */
