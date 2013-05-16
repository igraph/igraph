/* glpsql.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Author: Heinrich Schuchardt <xypron.glpk@gmx.de>.
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
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glpmpl.h"
#include "glpsql.h"

#ifdef ODBC_DLNAME
#define HAVE_ODBC
#define libodbc ODBC_DLNAME
#define h_odbc (get_env_ptr()->h_odbc)
#endif

#ifdef MYSQL_DLNAME
#define HAVE_MYSQL
#define libmysql MYSQL_DLNAME
#define h_mysql (get_env_ptr()->h_mysql)
#endif

static void *db_iodbc_open_int(TABDCA *dca, int mode, const char
      **sqllines);
static void *db_mysql_open_int(TABDCA *dca, int mode, const char
      **sqllines);

/**********************************************************************/

#if defined(HAVE_ODBC) || defined(HAVE_MYSQL)

#define SQL_FIELD_MAX 100
/* maximal field count */

#define SQL_FDLEN_MAX 255
/* maximal field length */

/***********************************************************************
*  NAME
*
*  args_concat - concatenate arguments
*
*  SYNOPSIS
*
*  static char **args_concat(TABDCA *dca);
*
*  DESCRIPTION
*
*  The arguments passed in dca are SQL statements. A SQL statement may
*  be split over multiple arguments. The last argument of a SQL
*  statement will be terminated with a semilocon. Each SQL statement is
*  merged into a single zero terminated string. Boundaries between
*  arguments are replaced by space.
*
*  RETURNS
*
*  Buffer with SQL statements */

static char **args_concat(TABDCA *dca)
{
   const char  *arg;
   int          i;
   int          j;
   int          j0;
   int          j1;
   int          len;
   int          lentot;
   int          narg;
   int          nline = 0;
   void        *ret;
   char       **sqllines = NULL;

   narg = mpl_tab_num_args(dca);
   /* The SQL statements start with argument 3. */
   if (narg < 3)
      return NULL;
   /* Count the SQL statements */
   for (j = 3; j <= narg; j++)
   {
      arg = mpl_tab_get_arg(dca, j);
      len = strlen(arg);
      if (arg[len-1] == ';' || j == narg)
        nline ++;
   }
   /* Allocate string buffer. */
   sqllines = (char **) xmalloc((nline+1) * sizeof(char **));
   /* Join arguments */
   sqllines[0] = NULL;
   j0     = 3;
   i      = 0;
   lentot = 0;
   for (j = 3; j <= narg; j++)
   {
      arg = mpl_tab_get_arg(dca, j);
      len = strlen(arg);
      lentot += len;
      if (arg[len-1] == ';' || j == narg)
      {  /* Join arguments for a single SQL statement */
         sqllines[i] = xmalloc(lentot+1);
         sqllines[i+1] = NULL;
         sqllines[i][0] = 0x00;
         for (j1 = j0; j1 <= j; j1++)
         {  if(j1>j0)
               strcat(sqllines[i], " ");
            strcat(sqllines[i], mpl_tab_get_arg(dca, j1));
         }
         len = strlen(sqllines[i]);
         if (sqllines[i][len-1] == ';')
            sqllines[i][len-1] = 0x00;
         j0 = j+1;
         i++;
         lentot = 0;
      }
   }
   return sqllines;
}

/***********************************************************************
*  NAME
*
*  free_buffer - free multiline string buffer
*
*  SYNOPSIS
*
*  static void free_buffer(char **buf);
*
*  DESCRIPTION
*
*  buf is a list of strings terminated by NULL.
*  The memory for the strings and for the list is released. */

static void free_buffer(char **buf)
{  int i;

   for(i = 0; buf[i] != NULL; i++)
      xfree(buf[i]);
   xfree(buf);
}

static int db_escaped_string_length(const char* from)
/* length of escaped string */
{
   int         count;
   const char *pointer;

    for (pointer = from, count = 0; *pointer != (char) '\0'; pointer++,
         count++)
    {
      switch (*pointer)
      {
         case '\'':
            count++;
            break;
      }
    }

    return count;
}

static int db_escape_string (char *to, const char *from)
/* escape string*/
{
   const char *source = from;
   char *target = to;
   unsigned int remaining;

   remaining = strlen(from);

   if (to == NULL)
     to = (char *) (from + remaining);

   while (remaining > 0)
   {
      switch (*source)
      {
         case '\'':
            *target = '\'';
            target++;
            *target = '\'';
            break;

         default:
            *target = *source;
            }
      source++;
      target++;
      remaining--;
      }

   /* Write the terminating NUL character. */
   *target = '\0';

   return target - to;
}

static char *db_generate_select_stmt(TABDCA *dca)
/* generate select statement */
{
   char        *arg;
   char const  *field;
   char        *query;
   int          j;
   int          narg;
   int          nf;
   int          total;

   total = 50;
   nf = mpl_tab_num_flds(dca);
   narg = mpl_tab_num_args(dca);
   for (j=1; j <= nf && j <= SQL_FIELD_MAX; j++)
   {
      field = mpl_tab_get_name(dca, j);
      total += strlen(field);
      total += 2;
   }
   arg = (char *) mpl_tab_get_arg(dca, narg);
   total += strlen(arg);
   query = xmalloc( total * sizeof(char));
   strcpy (query, "SELECT ");
   for (j=1; j <= nf && j <= SQL_FIELD_MAX; j++)
   {
      field = mpl_tab_get_name(dca, j);
      strcat(query, field);
      if ( j < nf )
         strcat(query, ", ");
   }
   strcat(query, " FROM ");
   strcat(query, arg);
   return query;
}

static char *db_generate_insert_stmt(TABDCA *dca)
/* generate insert statement */
{
   char        *arg;
   char const  *field;
   char        *query;
   int          j;
   int          narg;
   int          nf;
   int          total;

   total = 50;
   nf = mpl_tab_num_flds(dca);
   narg = mpl_tab_num_args(dca);
   for (j=1; j <= nf && j <= SQL_FIELD_MAX; j++)
   {
      field = mpl_tab_get_name(dca, j);
      total += strlen(field);
      total += 5;
   }
   arg = (char *) mpl_tab_get_arg(dca, narg);
   total += strlen(arg);
   query = xmalloc( (total+1) * sizeof(char));
   strcpy (query, "INSERT INTO ");
   strcat(query, arg);
   strcat(query, " ( ");
   for (j=1; j <= nf && j <= SQL_FIELD_MAX; j++)
   {
      field = mpl_tab_get_name(dca, j);
      strcat(query, field);
      if ( j < nf )
         strcat(query, ", ");
   }
   strcat(query, " ) VALUES ( ");
   for (j=1; j <= nf && j <= SQL_FIELD_MAX; j++)
   {
      strcat(query, "?");
      if ( j < nf )
         strcat(query, ", ");
   }
   strcat(query, " )");
   return query;
}

#endif

/**********************************************************************/

#ifndef HAVE_ODBC

void *db_iodbc_open(TABDCA *dca, int mode)
{     xassert(dca == dca);
      xassert(mode == mode);
      xprintf("iODBC table driver not supported\n");
      return NULL;
}

int db_iodbc_read(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

int db_iodbc_write(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

int db_iodbc_close(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

#else

#if defined(__CYGWIN__) || defined(__MINGW32__) || defined(__WOE__)
#include <windows.h>
#endif

#include <sql.h>
#include <sqlext.h>

struct db_odbc
{
   int              mode;         /*'R' = Read, 'W' = Write*/
   SQLHDBC          hdbc;         /*connection handle*/
   SQLHENV          henv;         /*environment handle*/
   SQLHSTMT         hstmt;        /*statement handle*/
   SQLSMALLINT      nresultcols;  /* columns in result*/
   SQLULEN          collen[SQL_FIELD_MAX+1];
   SQLLEN           outlen[SQL_FIELD_MAX+1];
   SQLSMALLINT      coltype[SQL_FIELD_MAX+1];
   SQLCHAR          data[SQL_FIELD_MAX+1][SQL_FDLEN_MAX+1];
   SQLCHAR          colname[SQL_FIELD_MAX+1][SQL_FDLEN_MAX+1];
   int              isnumeric[SQL_FIELD_MAX+1];
   int              nf;
   /* number of fields in the csv file */
   int              ref[1+SQL_FIELD_MAX];
   /* ref[k] = k', if k-th field of the csv file corresponds to
      k'-th field in the table statement; if ref[k] = 0, k-th field
      of the csv file is ignored */
   SQLCHAR         *query;
   /* query generated by db_iodbc_open */
};

SQLRETURN SQL_API dl_SQLAllocHandle (
   SQLSMALLINT           HandleType,
   SQLHANDLE             InputHandle,
   SQLHANDLE            *OutputHandle)
{
      typedef SQLRETURN SQL_API ep_SQLAllocHandle(
         SQLSMALLINT           HandleType,
         SQLHANDLE             InputHandle,
         SQLHANDLE            *OutputHandle);

      ep_SQLAllocHandle *fn;
      fn = (ep_SQLAllocHandle *) xdlsym(h_odbc, "SQLAllocHandle");
      xassert(fn != NULL);
      return (*fn)(HandleType, InputHandle, OutputHandle);
}

SQLRETURN SQL_API dl_SQLBindCol (
   SQLHSTMT              StatementHandle,
   SQLUSMALLINT          ColumnNumber,
   SQLSMALLINT           TargetType,
   SQLPOINTER            TargetValue,
   SQLLEN                BufferLength,
   SQLLEN               *StrLen_or_Ind)
{
      typedef SQLRETURN SQL_API ep_SQLBindCol(
         SQLHSTMT              StatementHandle,
         SQLUSMALLINT          ColumnNumber,
         SQLSMALLINT           TargetType,
         SQLPOINTER            TargetValue,
         SQLLEN                BufferLength,
         SQLLEN               *StrLen_or_Ind);
      ep_SQLBindCol *fn;
      fn = (ep_SQLBindCol *) xdlsym(h_odbc, "SQLBindCol");
      xassert(fn != NULL);
      return (*fn)(StatementHandle, ColumnNumber, TargetType,
         TargetValue, BufferLength, StrLen_or_Ind);
}

SQLRETURN SQL_API dl_SQLCloseCursor (
   SQLHSTMT              StatementHandle)
{
      typedef SQLRETURN SQL_API ep_SQLCloseCursor (
         SQLHSTMT              StatementHandle);

      ep_SQLCloseCursor *fn;
      fn = (ep_SQLCloseCursor *) xdlsym(h_odbc, "SQLCloseCursor");
      xassert(fn != NULL);
      return (*fn)(StatementHandle);
}


SQLRETURN SQL_API dl_SQLDisconnect (
   SQLHDBC               ConnectionHandle)
{
      typedef SQLRETURN SQL_API ep_SQLDisconnect(
         SQLHDBC               ConnectionHandle);

      ep_SQLDisconnect *fn;
      fn = (ep_SQLDisconnect *) xdlsym(h_odbc, "SQLDisconnect");
      xassert(fn != NULL);
      return (*fn)(ConnectionHandle);
}

SQLRETURN SQL_API dl_SQLDriverConnect (
   SQLHDBC               hdbc,
   SQLHWND               hwnd,
   SQLCHAR              *szConnStrIn,
   SQLSMALLINT           cbConnStrIn,
   SQLCHAR              *szConnStrOut,
   SQLSMALLINT           cbConnStrOutMax,
   SQLSMALLINT          *pcbConnStrOut,
   SQLUSMALLINT          fDriverCompletion)
{
      typedef SQLRETURN SQL_API ep_SQLDriverConnect(
         SQLHDBC               hdbc,
         SQLHWND               hwnd,
         SQLCHAR             * szConnStrIn,
         SQLSMALLINT           cbConnStrIn,
         SQLCHAR             * szConnStrOut,
         SQLSMALLINT           cbConnStrOutMax,
         SQLSMALLINT         * pcbConnStrOut,
         SQLUSMALLINT          fDriverCompletion);

      ep_SQLDriverConnect *fn;
      fn = (ep_SQLDriverConnect *) xdlsym(h_odbc, "SQLDriverConnect");
      xassert(fn != NULL);
      return (*fn)(hdbc, hwnd, szConnStrIn, cbConnStrIn, szConnStrOut,
         cbConnStrOutMax, pcbConnStrOut, fDriverCompletion);
}

SQLRETURN SQL_API dl_SQLEndTran (
   SQLSMALLINT           HandleType,
   SQLHANDLE             Handle,
   SQLSMALLINT           CompletionType)
{
      typedef SQLRETURN SQL_API ep_SQLEndTran (
         SQLSMALLINT           HandleType,
         SQLHANDLE             Handle,
         SQLSMALLINT           CompletionType);

      ep_SQLEndTran *fn;
      fn = (ep_SQLEndTran *) xdlsym(h_odbc, "SQLEndTran");
      xassert(fn != NULL);
      return (*fn)(HandleType, Handle, CompletionType);
}

SQLRETURN SQL_API dl_SQLExecDirect (
   SQLHSTMT              StatementHandle,
   SQLCHAR             * StatementText,
   SQLINTEGER            TextLength)
{
      typedef SQLRETURN SQL_API ep_SQLExecDirect (
         SQLHSTMT              StatementHandle,
         SQLCHAR             * StatementText,
         SQLINTEGER            TextLength);

      ep_SQLExecDirect *fn;
      fn = (ep_SQLExecDirect *) xdlsym(h_odbc, "SQLExecDirect");
      xassert(fn != NULL);
      return (*fn)(StatementHandle, StatementText, TextLength);
}

SQLRETURN SQL_API dl_SQLFetch (
   SQLHSTMT              StatementHandle)
{
      typedef SQLRETURN SQL_API ep_SQLFetch (
         SQLHSTMT              StatementHandle);

      ep_SQLFetch *fn;
      fn = (ep_SQLFetch*) xdlsym(h_odbc, "SQLFetch");
      xassert(fn != NULL);
      return (*fn)(StatementHandle);
}

SQLRETURN SQL_API dl_SQLFreeHandle (
   SQLSMALLINT           HandleType,
   SQLHANDLE             Handle)
{
      typedef SQLRETURN SQL_API ep_SQLFreeHandle (
         SQLSMALLINT           HandleType,
         SQLHANDLE             Handle);

      ep_SQLFreeHandle *fn;
      fn = (ep_SQLFreeHandle *) xdlsym(h_odbc, "SQLFreeHandle");
      xassert(fn != NULL);
      return (*fn)(HandleType, Handle);
}

SQLRETURN SQL_API dl_SQLDescribeCol (
   SQLHSTMT              StatementHandle,
   SQLUSMALLINT          ColumnNumber,
   SQLCHAR             * ColumnName,
   SQLSMALLINT           BufferLength,
   SQLSMALLINT         * NameLength,
   SQLSMALLINT         * DataType,
   SQLULEN             * ColumnSize,
   SQLSMALLINT         * DecimalDigits,
   SQLSMALLINT         * Nullable)
{
      typedef SQLRETURN SQL_API ep_SQLDescribeCol (
         SQLHSTMT              StatementHandle,
         SQLUSMALLINT          ColumnNumber,
         SQLCHAR              *ColumnName,
         SQLSMALLINT           BufferLength,
         SQLSMALLINT          *NameLength,
         SQLSMALLINT          *DataType,
         SQLULEN              *ColumnSize,
         SQLSMALLINT          *DecimalDigits,
         SQLSMALLINT          *Nullable);

      ep_SQLDescribeCol *fn;
      fn = (ep_SQLDescribeCol *) xdlsym(h_odbc, "SQLDescribeCol");
      xassert(fn != NULL);
      return (*fn)(StatementHandle, ColumnNumber, ColumnName,
         BufferLength, NameLength,
         DataType, ColumnSize, DecimalDigits, Nullable);
}

SQLRETURN SQL_API dl_SQLGetDiagRec (
   SQLSMALLINT           HandleType,
   SQLHANDLE             Handle,
   SQLSMALLINT           RecNumber,
   SQLCHAR              *Sqlstate,
   SQLINTEGER           *NativeError,
   SQLCHAR              *MessageText,
   SQLSMALLINT           BufferLength,
   SQLSMALLINT          *TextLength)
{
      typedef SQLRETURN SQL_API ep_SQLGetDiagRec (
         SQLSMALLINT           HandleType,
         SQLHANDLE             Handle,
         SQLSMALLINT           RecNumber,
         SQLCHAR              *Sqlstate,
         SQLINTEGER           *NativeError,
         SQLCHAR              *MessageText,
         SQLSMALLINT           BufferLength,
         SQLSMALLINT          *TextLength);

      ep_SQLGetDiagRec *fn;
      fn = (ep_SQLGetDiagRec *) xdlsym(h_odbc, "SQLGetDiagRec");
      xassert(fn != NULL);
      return (*fn)(HandleType, Handle, RecNumber, Sqlstate,
         NativeError, MessageText, BufferLength, TextLength);
}

SQLRETURN SQL_API dl_SQLGetInfo (
   SQLHDBC               ConnectionHandle,
   SQLUSMALLINT          InfoType,
   SQLPOINTER            InfoValue,
   SQLSMALLINT           BufferLength,
   SQLSMALLINT          *StringLength)
{
      typedef SQLRETURN SQL_API ep_SQLGetInfo (
         SQLHDBC               ConnectionHandle,
         SQLUSMALLINT          InfoType,
         SQLPOINTER            InfoValue,
         SQLSMALLINT           BufferLength,
         SQLSMALLINT          *StringLength);

      ep_SQLGetInfo *fn;
      fn = (ep_SQLGetInfo *) xdlsym(h_odbc, "SQLGetInfo");
      xassert(fn != NULL);
      return (*fn)(ConnectionHandle, InfoType, InfoValue, BufferLength,
         StringLength);
}

SQLRETURN SQL_API dl_SQLNumResultCols (
   SQLHSTMT              StatementHandle,
   SQLSMALLINT          *ColumnCount)
{
      typedef SQLRETURN SQL_API ep_SQLNumResultCols (
         SQLHSTMT              StatementHandle,
         SQLSMALLINT          *ColumnCount);

      ep_SQLNumResultCols *fn;
      fn = (ep_SQLNumResultCols *) xdlsym(h_odbc, "SQLNumResultCols");
      xassert(fn != NULL);
      return (*fn)(StatementHandle, ColumnCount);
}

SQLRETURN SQL_API dl_SQLSetConnectAttr (
   SQLHDBC               ConnectionHandle,
   SQLINTEGER            Attribute,
   SQLPOINTER            Value,
   SQLINTEGER            StringLength)
{
      typedef SQLRETURN SQL_API ep_SQLSetConnectAttr (
         SQLHDBC               ConnectionHandle,
         SQLINTEGER            Attribute,
         SQLPOINTER            Value,
         SQLINTEGER            StringLength);

      ep_SQLSetConnectAttr *fn;
     fn = (ep_SQLSetConnectAttr *) xdlsym(h_odbc, "SQLSetConnectAttr");
      xassert(fn != NULL);
      return (*fn)(ConnectionHandle, Attribute, Value, StringLength);
}

SQLRETURN SQL_API dl_SQLSetEnvAttr (
   SQLHENV               EnvironmentHandle,
   SQLINTEGER            Attribute,
   SQLPOINTER            Value,
   SQLINTEGER            StringLength)
{
      typedef SQLRETURN SQL_API ep_SQLSetEnvAttr (
         SQLHENV               EnvironmentHandle,
         SQLINTEGER            Attribute,
         SQLPOINTER            Value,
         SQLINTEGER            StringLength);

      ep_SQLSetEnvAttr *fn;
      fn = (ep_SQLSetEnvAttr *) xdlsym(h_odbc, "SQLSetEnvAttr");
      xassert(fn != NULL);
      return (*fn)(EnvironmentHandle, Attribute, Value, StringLength);
}

static void extract_error(
   char *fn,
   SQLHANDLE handle,
   SQLSMALLINT type);

static int is_numeric(
    SQLSMALLINT coltype);

/***********************************************************************
*  NAME
*
*  db_iodbc_open - open connection to ODBC data base
*
*  SYNOPSIS
*
*  #include "glpsql.h"
*  void *db_iodbc_open(TABDCA *dca, int mode);
*
*  DESCRIPTION
*
*  The routine db_iodbc_open opens a connection to an ODBC data base.
*  It then executes the sql statements passed.
*
*  In the case of table read the SELECT statement is executed.
*
*  In the case of table write the INSERT statement is prepared.
*  RETURNS
*
*  The routine returns a pointer to data storage area created. */
void *db_iodbc_open(TABDCA *dca, int mode)
{  void  *ret;
   char **sqllines;

   sqllines = args_concat(dca);
   if (sqllines == NULL)
   {  xprintf("Missing arguments in table statement.\n"
              "Please, supply table driver, dsn, and query.\n");
      return NULL;
   }
   ret = db_iodbc_open_int(dca, mode, (const char **) sqllines);
   free_buffer(sqllines);
   return ret;
}

static void *db_iodbc_open_int(TABDCA *dca, int mode, const char
   **sqllines)
{
   struct db_odbc    *sql;
   SQLRETURN          ret;
   SQLCHAR FAR       *dsn;
   SQLCHAR            info[256];
   SQLSMALLINT        colnamelen;
   SQLSMALLINT        nullable;
   SQLSMALLINT        scale;
   const char        *arg;
   int                narg;
   int                i, j;
   int                total;

   if (libodbc == NULL)
   {
      xprintf("No loader for shared ODBC library available\n");
      return NULL;
   }

   if (h_odbc == NULL)
   {
      h_odbc = xdlopen(libodbc);
      if (h_odbc == NULL)
      {  xprintf("unable to open library %s\n", libodbc);
         xprintf("%s\n", xerrmsg());
         return NULL;
      }
   }

   sql = (struct db_odbc *) xmalloc(sizeof(struct db_odbc));
   if (sql == NULL)
         return NULL;

   sql->mode  = mode;
   sql->hdbc  = NULL;
   sql->henv  = NULL;
   sql->hstmt = NULL;
   sql->query = NULL;
   narg = mpl_tab_num_args(dca);

   dsn = (SQLCHAR FAR *) mpl_tab_get_arg(dca, 2);
   /* allocate an environment handle */
   ret = dl_SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE,
      &(sql->henv));
   /* set attribute to enable application to run as ODBC 3.0
      application */
   ret = dl_SQLSetEnvAttr(sql->henv, SQL_ATTR_ODBC_VERSION,
      (void *) SQL_OV_ODBC3, 0);
   /* allocate a connection handle */
   ret = dl_SQLAllocHandle(SQL_HANDLE_DBC, sql->henv, &(sql->hdbc));
   /* connect */
   ret = dl_SQLDriverConnect(sql->hdbc, NULL, dsn, SQL_NTS, NULL, 0,
      NULL, SQL_DRIVER_COMPLETE);
   if (SQL_SUCCEEDED(ret))
   {  /* output information about data base connection */
      xprintf("Connected to ");
      dl_SQLGetInfo(sql->hdbc, SQL_DBMS_NAME, (SQLPOINTER)info,
         sizeof(info), NULL);
      xprintf("%s ", info);
      dl_SQLGetInfo(sql->hdbc, SQL_DBMS_VER, (SQLPOINTER)info,
         sizeof(info), NULL);
      xprintf("%s - ", info);
      dl_SQLGetInfo(sql->hdbc, SQL_DATABASE_NAME, (SQLPOINTER)info,
         sizeof(info), NULL);
      xprintf("%s\n", info);
   }
   else
   {  /* describe error */
      xprintf("Failed to connect\n");
      extract_error("SQLDriverConnect", sql->hdbc, SQL_HANDLE_DBC);
      dl_SQLFreeHandle(SQL_HANDLE_DBC, sql->hdbc);
      dl_SQLFreeHandle(SQL_HANDLE_ENV, sql->henv);
      xfree(sql);
      return NULL;
   }
   /* set AUTOCOMMIT on*/
   ret = dl_SQLSetConnectAttr(sql->hdbc, SQL_ATTR_AUTOCOMMIT,
      (SQLPOINTER)SQL_AUTOCOMMIT_ON, 0);
   /* allocate a statement handle */
   ret = dl_SQLAllocHandle(SQL_HANDLE_STMT, sql->hdbc, &(sql->hstmt));

   /* initialization queries */
   for(j = 0; sqllines[j+1] != NULL; j++)
   {
      sql->query = (SQLCHAR *) sqllines[j];
      xprintf("%s\n", sql->query);
      ret = dl_SQLExecDirect(sql->hstmt, sql->query, SQL_NTS);
      switch (ret)
      {
         case SQL_SUCCESS:
         case SQL_SUCCESS_WITH_INFO:
         case SQL_NO_DATA_FOUND:
            break;
         default:
            xprintf("db_iodbc_open: Query\n\"%s\"\nfailed.\n",
               sql->query);
            extract_error("SQLExecDirect", sql->hstmt, SQL_HANDLE_STMT);
            dl_SQLFreeHandle(SQL_HANDLE_STMT, sql->hstmt);
            dl_SQLDisconnect(sql->hdbc);
            dl_SQLFreeHandle(SQL_HANDLE_DBC, sql->hdbc);
            dl_SQLFreeHandle(SQL_HANDLE_ENV, sql->henv);
            xfree(sql);
            return NULL;
      }
      /* commit statement */
      dl_SQLEndTran(SQL_HANDLE_ENV, sql->henv, SQL_COMMIT);
   }

   if ( sql->mode == 'R' )
   {  sql->nf = mpl_tab_num_flds(dca);
      for(j = 0; sqllines[j] != NULL; j++)
         arg = sqllines[j];
      total = strlen(arg);
      if (total > 7 && 0 == strncmp(arg, "SELECT ", 7))
      {
         total = strlen(arg);
         sql->query = xmalloc( (total+1) * sizeof(char));
         strcpy (sql->query, arg);
      }
      else
      {
         sql->query = db_generate_select_stmt(dca);
      }
      xprintf("%s\n", sql->query);
      if (dl_SQLExecDirect(sql->hstmt, sql->query, SQL_NTS) !=
         SQL_SUCCESS)
      {
         xprintf("db_iodbc_open: Query\n\"%s\"\nfailed.\n", sql->query);
         extract_error("SQLExecDirect", sql->hstmt, SQL_HANDLE_STMT);
         dl_SQLFreeHandle(SQL_HANDLE_STMT, sql->hstmt);
         dl_SQLDisconnect(sql->hdbc);
         dl_SQLFreeHandle(SQL_HANDLE_DBC, sql->hdbc);
         dl_SQLFreeHandle(SQL_HANDLE_ENV, sql->henv);
         xfree(sql->query);
            xfree(sql);
         return NULL;
      }
      xfree(sql->query);
      /* determine number of result columns */
      ret = dl_SQLNumResultCols(sql->hstmt, &sql->nresultcols);
      total = sql->nresultcols;
      if (total > SQL_FIELD_MAX)
      {  xprintf("db_iodbc_open: Too many fields (> %d) in query.\n"
            "\"%s\"\n", SQL_FIELD_MAX, sql->query);
         dl_SQLFreeHandle(SQL_HANDLE_STMT, sql->hstmt);
         dl_SQLDisconnect(sql->hdbc);
         dl_SQLFreeHandle(SQL_HANDLE_DBC, sql->hdbc);
         dl_SQLFreeHandle(SQL_HANDLE_ENV, sql->henv);
         xfree(sql->query);
         return NULL;
      }
      for (i = 1; i <= total; i++)
      {  /* return a set of attributes for a column */
         ret = dl_SQLDescribeCol(sql->hstmt, (SQLSMALLINT) i,
            sql->colname[i], SQL_FDLEN_MAX,
            &colnamelen, &(sql->coltype[i]), &(sql->collen[i]), &scale,
            &nullable);
         sql->isnumeric[i] = is_numeric(sql->coltype[i]);
         /* bind columns to program vars, converting all types to CHAR*/
         dl_SQLBindCol(sql->hstmt, i, SQL_CHAR, sql->data[i],
            SQL_FDLEN_MAX, &(sql->outlen[i]));
         for (j = sql->nf; j >= 1; j--)
         {  if (strcmp(mpl_tab_get_name(dca, j), sql->colname[i]) == 0)
            break;
         }
         sql->ref[i] = j;
      }
   }
   else if ( sql->mode == 'W' )
   {  for(j = 0; sqllines[j] != NULL; j++)
         arg = sqllines[j];
      if (  NULL != strchr(arg, '?') )
      {
         total = strlen(arg);
         sql->query = xmalloc( (total+1) * sizeof(char));
         strcpy (sql->query, arg);
         }
      else
      {
         sql->query = db_generate_insert_stmt(dca);
      }
      xprintf("%s\n", sql->query);
   }
   return sql;
}

int db_iodbc_read(TABDCA *dca, void *link)
{
   struct db_odbc  *sql;
   SQLRETURN        ret;
   char             buf[SQL_FDLEN_MAX+1];
   int              i;
   int              len;
   double           num;

   sql = (struct db_odbc *) link;

   xassert(sql != NULL);
   xassert(sql->mode == 'R');

   ret=dl_SQLFetch(sql->hstmt);
   if (ret== SQL_ERROR)
      return -1;
   if (ret== SQL_NO_DATA_FOUND)
      return -1; /*EOF*/
   for (i=1; i <= sql->nresultcols; i++)
   {
      if (sql->ref[i] > 0)
      {
         len = sql->outlen[i];
         if (len != SQL_NULL_DATA)
         {
            if (len > SQL_FDLEN_MAX)
               len = SQL_FDLEN_MAX;
            else if (len < 0)
               len = 0;
            strncpy(buf, (const char *) sql->data[i], len);
            buf[len] = 0x00;
            if (0 != (sql->isnumeric[i]))
            {  strspx(buf); /* remove spaces*/
               if (str2num(buf, &num) != 0)
               {  xprintf("'%s' cannot be converted to a number.\n",
                     buf);
                  return 1;
               }
               mpl_tab_set_num(dca, sql->ref[i], num);
            }
            else
            {  mpl_tab_set_str(dca, sql->ref[i], strtrim(buf));
            }
         }
      }
   }
   return 0;
}

int db_iodbc_write(TABDCA *dca, void *link)
{
   struct db_odbc  *sql;
   char            *part;
   char            *query;
   char            *template;
   char             num[50];
   int              k;
   int              len;
   int              nf;

   sql = (struct db_odbc *) link;
   xassert(sql != NULL);
   xassert(sql->mode == 'W');

   len      = strlen(sql->query);
   template = (char *) xmalloc( (len + 1) * sizeof(char) );
   strcpy(template, sql->query);

   nf = mpl_tab_num_flds(dca);
   for (k = 1; k <= nf; k++)
   {     switch (mpl_tab_get_type(dca, k))
      {  case 'N':
            len += 20;
            break;
         case 'S':
            len += db_escaped_string_length(mpl_tab_get_str(dca, k));
            len += 2;
            break;
              default:
                        xassert(dca != dca);
         }
   }
   query = xmalloc( (len + 1 ) * sizeof(char) );
   query[0] = 0x00;
   for (k = 1, part = strtok (template, "?"); (part != NULL);
      part = strtok (NULL, "?"), k++)
   {
      if (k > nf) break;
      strcat( query, part );
      switch (mpl_tab_get_type(dca, k))
      {  case 'N':
#if 0 /* 02/XI-2010 by xypron */
            sprintf(num, "%-18g",mpl_tab_get_num(dca, k));
#else
            sprintf(num, "%.*g", DBL_DIG, mpl_tab_get_num(dca, k));
#endif
            strcat( query, num );
            break;
         case 'S':
            strcat( query, "'");
            db_escape_string( query + strlen(query),
               mpl_tab_get_str(dca, k) );
            strcat( query, "'");
            break;
              default:
                        xassert(dca != dca);
         }
   }
   if (part != NULL)
      strcat(query, part);
   if (dl_SQLExecDirect(sql->hstmt, (SQLCHAR *) query, SQL_NTS)
      != SQL_SUCCESS)
   {
      xprintf("db_iodbc_write: Query\n\"%s\"\nfailed.\n", query);
      extract_error("SQLExecDirect", sql->hdbc, SQL_HANDLE_DBC);
      xfree(query);
      xfree(template);
      return 1;
      }

   xfree(query);
   xfree(template);
   return 0;
}

int db_iodbc_close(TABDCA *dca, void *link)
{
   struct db_odbc *sql;

   sql = (struct db_odbc *) link;
   xassert(sql != NULL);
   /* Commit */
   if ( sql->mode == 'W' )
      dl_SQLEndTran(SQL_HANDLE_ENV, sql->henv, SQL_COMMIT);
   if ( sql->mode == 'R' )
      dl_SQLCloseCursor(sql->hstmt);

   dl_SQLFreeHandle(SQL_HANDLE_STMT, sql->hstmt);
   dl_SQLDisconnect(sql->hdbc);
   dl_SQLFreeHandle(SQL_HANDLE_DBC, sql->hdbc);
   dl_SQLFreeHandle(SQL_HANDLE_ENV, sql->henv);
   if ( sql->mode == 'W' )
      xfree(sql->query);
   xfree(sql);
   dca->link = NULL;
   return 0;
}

static void extract_error(
   char *fn,
   SQLHANDLE handle,
   SQLSMALLINT type)
{
   SQLINTEGER   i = 0;
   SQLINTEGER   native;
   SQLCHAR   state[ 7 ];
   SQLCHAR   text[256];
   SQLSMALLINT  len;
   SQLRETURN    ret;

   xprintf("\nThe driver reported the following diagnostics whilst "
      "running %s\n", fn);

   do
   {
      ret = dl_SQLGetDiagRec(type, handle, ++i, state, &native, text,
         sizeof(text), &len );
      if (SQL_SUCCEEDED(ret))
         xprintf("%s:%ld:%ld:%s\n", state, i, native, text);
   }
   while( ret == SQL_SUCCESS );
}

static int is_numeric(SQLSMALLINT coltype)
{
   int ret = 0;
   switch (coltype)
   {
      case SQL_DECIMAL:
      case SQL_NUMERIC:
      case SQL_SMALLINT:
      case SQL_INTEGER:
      case SQL_REAL:
      case SQL_FLOAT:
      case SQL_DOUBLE:
      case SQL_TINYINT:
      case SQL_BIGINT:
         ret = 1;
         break;
   }
   return ret;
}

#endif

/**********************************************************************/

#ifndef HAVE_MYSQL

void *db_mysql_open(TABDCA *dca, int mode)
{     xassert(dca == dca);
      xassert(mode == mode);
      xprintf("MySQL table driver not supported\n");
      return NULL;
}

int db_mysql_read(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

int db_mysql_write(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

int db_mysql_close(TABDCA *dca, void *link)
{     xassert(dca != dca);
      xassert(link != link);
      return 0;
}

#else

#if defined(__CYGWIN__) || defined(__MINGW32__) || defined(__WOE__)
#include <windows.h>
#endif

#ifdef __CYGWIN__
#define byte_defined 1
#endif

#include <my_global.h>
#include <my_sys.h>
#include <mysql.h>

struct db_mysql
{
   int              mode;  /*'R' = Read, 'W' = Write*/
   MYSQL           *con;   /*connection*/
   MYSQL_RES       *res;    /*result*/
   int              nf;
   /* number of fields in the csv file */
   int              ref[1+SQL_FIELD_MAX];
   /* ref[k] = k', if k-th field of the csv file corresponds to
      k'-th field in the table statement; if ref[k] = 0, k-th field
      of the csv file is ignored */
   char            *query;
   /* query generated by db_mysql_open */
};

void STDCALL dl_mysql_close(MYSQL *sock)
{
      typedef void STDCALL ep_mysql_close(MYSQL *sock);

      ep_mysql_close *fn;
      fn = (ep_mysql_close *) xdlsym(h_mysql, "mysql_close");
      xassert(fn != NULL);
      return (*fn)(sock);
}

const char * STDCALL dl_mysql_error(MYSQL *mysql)
{
      typedef const char * STDCALL ep_mysql_error(MYSQL *mysql);

      ep_mysql_error *fn;
      fn = (ep_mysql_error *) xdlsym(h_mysql, "mysql_error");
      xassert(fn != NULL);
      return (*fn)(mysql);
}

MYSQL_FIELD * STDCALL dl_mysql_fetch_fields(MYSQL_RES *res)
{
      typedef MYSQL_FIELD * STDCALL
         ep_mysql_fetch_fields(MYSQL_RES *res);

      ep_mysql_fetch_fields *fn;
   fn = (ep_mysql_fetch_fields *) xdlsym(h_mysql, "mysql_fetch_fields");
      xassert(fn != NULL);
      return (*fn)(res);
}

unsigned long * STDCALL dl_mysql_fetch_lengths(MYSQL_RES *result)
{
      typedef unsigned long * STDCALL
         ep_mysql_fetch_lengths(MYSQL_RES *result);

      ep_mysql_fetch_lengths *fn;
      fn = (ep_mysql_fetch_lengths *) xdlsym(h_mysql,
         "mysql_fetch_lengths");
      xassert(fn != NULL);
      return (*fn)(result);
}

MYSQL_ROW STDCALL dl_mysql_fetch_row(MYSQL_RES *result)
{
      typedef MYSQL_ROW STDCALL ep_mysql_fetch_row(MYSQL_RES *result);

      ep_mysql_fetch_row *fn;
      fn = (ep_mysql_fetch_row *) xdlsym(h_mysql, "mysql_fetch_row");
      xassert(fn != NULL);
      return (*fn)(result);
}

unsigned int STDCALL dl_mysql_field_count(MYSQL *mysql)
{
      typedef unsigned int STDCALL ep_mysql_field_count(MYSQL *mysql);

      ep_mysql_field_count *fn;
     fn = (ep_mysql_field_count *) xdlsym(h_mysql, "mysql_field_count");
      xassert(fn != NULL);
      return (*fn)(mysql);
}

MYSQL * STDCALL dl_mysql_init(MYSQL *mysql)
{
      typedef MYSQL * STDCALL ep_mysql_init(MYSQL *mysql);

      ep_mysql_init *fn;
      fn = (ep_mysql_init *) xdlsym(h_mysql, "mysql_init");
      xassert(fn != NULL);
      return (*fn)(mysql);
}

unsigned int STDCALL dl_mysql_num_fields(MYSQL_RES *res)
{
      typedef unsigned int STDCALL ep_mysql_num_fields(MYSQL_RES *res);

      ep_mysql_num_fields *fn;
      fn = (ep_mysql_num_fields *) xdlsym(h_mysql, "mysql_num_fields");
      xassert(fn != NULL);
      return (*fn)(res);
}

int STDCALL dl_mysql_query(MYSQL *mysql, const char *q)
{
      typedef int STDCALL ep_mysql_query(MYSQL *mysql, const char *q);

      ep_mysql_query *fn;
      fn = (ep_mysql_query *) xdlsym(h_mysql, "mysql_query");
      xassert(fn != NULL);
      return (*fn)(mysql, q);
}

MYSQL * STDCALL dl_mysql_real_connect(MYSQL *mysql, const char *host,
                                           const char *user,
                                           const char *passwd,
                                           const char *db,
                                           unsigned int port,
                                           const char *unix_socket,
                                           unsigned long clientflag)
{
      typedef MYSQL * STDCALL ep_mysql_real_connect(MYSQL *mysql,
            const char *host,
            const char *user,
            const char *passwd,
            const char *db,
            unsigned int port,
            const char *unix_socket,
            unsigned long clientflag);

      ep_mysql_real_connect *fn;
      fn = (ep_mysql_real_connect *) xdlsym(h_mysql,
         "mysql_real_connect");
      xassert(fn != NULL);
      return (*fn)(mysql, host, user, passwd, db, port, unix_socket,
         clientflag);
}

MYSQL_RES * STDCALL dl_mysql_use_result(MYSQL *mysql)
{
      typedef MYSQL_RES * STDCALL ep_mysql_use_result(MYSQL *mysql);
      ep_mysql_use_result *fn;
      fn = (ep_mysql_use_result *) xdlsym(h_mysql, "mysql_use_result");
      xassert(fn != NULL);
      return (*fn)(mysql);
}

/***********************************************************************
*  NAME
*
*  db_mysql_open - open connection to ODBC data base
*
*  SYNOPSIS
*
*  #include "glpsql.h"
*  void *db_mysql_open(TABDCA *dca, int mode);
*
*  DESCRIPTION
*
*  The routine db_mysql_open opens a connection to a MySQL data base.
*  It then executes the sql statements passed.
*
*  In the case of table read the SELECT statement is executed.
*
*  In the case of table write the INSERT statement is prepared.
*  RETURNS
*
*  The routine returns a pointer to data storage area created. */

void *db_mysql_open(TABDCA *dca, int mode)
{  void  *ret;
   char **sqllines;

   sqllines = args_concat(dca);
   if (sqllines == NULL)
   {  xprintf("Missing arguments in table statement.\n"
              "Please, supply table driver, dsn, and query.\n");
      return NULL;
   }
   ret = db_mysql_open_int(dca, mode, (const char **) sqllines);
   free_buffer(sqllines);
   return ret;
}

static void *db_mysql_open_int(TABDCA *dca, int mode, const char
   **sqllines)
{
   struct db_mysql *sql = NULL;
   char            *arg = NULL;
   const char      *field;
   MYSQL_FIELD     *fields;
   char            *keyword;
   char            *value;
   char            *query;
   char            *dsn;
/* "Server=[server_name];Database=[database_name];UID=[username];*/
/* PWD=[password];Port=[port]"*/
   char            *server   = NULL;        /* Server */
   char            *user     = NULL;        /* UID */
   char            *password = NULL;        /* PWD */
   char            *database = NULL;        /* Database */
   unsigned int     port = 0;               /* Port */
   int              narg;
   int              i, j, total;

   if (libmysql == NULL)
   {
      xprintf("No loader for shared MySQL library available\n");
      return NULL;
   }

   if (h_mysql == NULL)
   {
      h_mysql = xdlopen(libmysql);
      if (h_mysql == NULL)
      {  xprintf("unable to open library %s\n", libmysql);
         xprintf("%s\n", xerrmsg());
         return NULL;
      }
   }

   sql = (struct db_mysql *) xmalloc(sizeof(struct db_mysql));
   if (sql == NULL)
         return NULL;
   sql->mode = mode;
   sql->res = NULL;
   sql->query = NULL;
   sql->nf = mpl_tab_num_flds(dca);

   narg = mpl_tab_num_args(dca);
   if (narg < 3 )
      xprintf("MySQL driver: string list too short \n");

   /* get connection string*/
   dsn = (char *) mpl_tab_get_arg(dca, 2);
      /* copy connection string*/
   i = strlen(dsn);
   i++;
   arg = xmalloc(i * sizeof(char));
   strcpy(arg, dsn);
   /*tokenize connection string*/
   for (i = 1, keyword = strtok (arg, "="); (keyword != NULL);
      keyword = strtok (NULL, "="), i++)
   {
         value = strtok (NULL, ";");
      if (value==NULL)
         {
            xprintf("db_mysql_open: Missing value for keyword %s\n",
               keyword);
            xfree(arg);
            xfree(sql);
            return NULL;
      }
      if (0 == strcmp(keyword, "Server"))
            server = value;
      else if (0 == strcmp(keyword, "Database"))
             database = value;
      else if (0 == strcmp(keyword, "UID"))
             user = value;
      else if (0 == strcmp(keyword, "PWD"))
             password = value;
      else if (0 == strcmp(keyword, "Port"))
             port = (unsigned int) atol(value);
   }
   /* Connect to database */
   sql->con = dl_mysql_init(NULL);
  if (!dl_mysql_real_connect(sql->con, server, user, password, database,
      port, NULL, 0))
   {
      xprintf("db_mysql_open: Connect failed\n");
      xprintf("%s\n", dl_mysql_error(sql->con));
      xfree(arg);
      xfree(sql);
      return NULL;
   }
   xfree(arg);

   for(j = 0; sqllines[j+1] != NULL; j++)
   {  query = (char *) sqllines[j];
      xprintf("%s\n", query);
      if (dl_mysql_query(sql->con, query))
      {
         xprintf("db_mysql_open: Query\n\"%s\"\nfailed.\n", query);
         xprintf("%s\n",dl_mysql_error(sql->con));
         dl_mysql_close(sql->con);
         xfree(sql);
         return NULL;
      }
   }

   if ( sql->mode == 'R' )
   {  sql->nf = mpl_tab_num_flds(dca);
      for(j = 0; sqllines[j] != NULL; j++)
         arg = (char *) sqllines[j];
      total = strlen(arg);
      if (total > 7 && 0 == strncmp(arg, "SELECT ", 7))
      {
         total = strlen(arg);
         query = xmalloc( (total+1) * sizeof(char));
         strcpy (query, arg);
      }
      else
      {
         query = db_generate_select_stmt(dca);
      }
      xprintf("%s\n", query);
      if (dl_mysql_query(sql->con, query))
      {
         xprintf("db_mysql_open: Query\n\"%s\"\nfailed.\n", query);
         xprintf("%s\n",dl_mysql_error(sql->con));
         dl_mysql_close(sql->con);
         xfree(query);
         xfree(sql);
         return NULL;
      }
      xfree(query);
      sql->res = dl_mysql_use_result(sql->con);
      if (sql->res)
      {
         /* create references between query results and table fields*/
         total = dl_mysql_num_fields(sql->res);
         if (total > SQL_FIELD_MAX)
         {  xprintf("db_mysql_open: Too many fields (> %d) in query.\n"
               "\"%s\"\n", SQL_FIELD_MAX, query);
            xprintf("%s\n",dl_mysql_error(sql->con));
            dl_mysql_close(sql->con);
            xfree(query);
                 xfree(sql);
            return NULL;
         }
         fields = dl_mysql_fetch_fields(sql->res);
         for (i = 1; i <= total; i++)
         {
               for (j = sql->nf; j >= 1; j--)
            {
               if (strcmp(mpl_tab_get_name(dca, j), fields[i-1].name)
                  == 0)
               break;
            }
            sql->ref[i] = j;
         }
      }
      else
      {
         if(dl_mysql_field_count(sql->con) == 0)
            {
            xprintf("db_mysql_open: Query was not a SELECT\n\"%s\"\n",
               query);
            xprintf("%s\n",dl_mysql_error(sql->con));
            xfree(query);
            xfree(sql);
            return NULL;
         }
         else
         {
            xprintf("db_mysql_open: Query\n\"%s\"\nfailed.\n", query);
            xprintf("%s\n",dl_mysql_error(sql->con));
            xfree(query);
            xfree(sql);
            return NULL;
         }
      }
   }
   else if ( sql->mode == 'W' )
   {  for(j = 0; sqllines[j] != NULL; j++)
         arg = (char *) sqllines[j];
      if (  NULL != strchr(arg, '?') )
      {
         total = strlen(arg);
         query = xmalloc( (total+1) * sizeof(char));
         strcpy (query, arg);
         }
      else
         query = db_generate_insert_stmt(dca);
      sql->query = query;
      xprintf("%s\n", query);
   }
   return sql;
}

int db_mysql_read(TABDCA *dca, void *link)
{  struct db_mysql *sql;
   char            buf[255+1];
   char            **row;
   unsigned long   *lengths;
   MYSQL_FIELD     *fields;
   double          num;
   int             len;
   unsigned long   num_fields;
   int             i;

   sql = (struct db_mysql *) link;

   xassert(sql != NULL);
   xassert(sql->mode == 'R');
   if (NULL == sql->res)
   {
      xprintf("db_mysql_read: no result set available");
      return 1;
   }
   if (NULL==(row = (char **)dl_mysql_fetch_row(sql->res))) {
       return -1; /*EOF*/
   }
   lengths = dl_mysql_fetch_lengths(sql->res);
   fields = dl_mysql_fetch_fields(sql->res);
   num_fields = dl_mysql_num_fields(sql->res);
   for (i=1; i <= num_fields; i++)
   {
      if (row[i-1] != NULL)
      {  len = (size_t) lengths[i-1];
         if (len > 255)
            len = 255;
         strncpy(buf, (const char *) row[i-1], len);
         buf[len] = 0x00;
         if (0 != (fields[i-1].flags & NUM_FLAG))
         {  strspx(buf); /* remove spaces*/
            if (str2num(buf, &num) != 0)
            {  xprintf("'%s' cannot be converted to a number.\n", buf);
               return 1;
            }
            if (sql->ref[i] > 0)
               mpl_tab_set_num(dca, sql->ref[i], num);
         }
         else
         {  if (sql->ref[i] > 0)
               mpl_tab_set_str(dca, sql->ref[i], strtrim(buf));
         }
      }
   }
   return 0;
}

int db_mysql_write(TABDCA *dca, void *link)
{
   struct db_mysql *sql;
   char            *part;
   char            *query;
   char            *template;
   char             num[50];
   int              k;
   int              len;
   int              nf;

   sql = (struct db_mysql *) link;
   xassert(sql != NULL);
   xassert(sql->mode == 'W');

   len      = strlen(sql->query);
   template = (char *) xmalloc( (len + 1) * sizeof(char) );
   strcpy(template, sql->query);

   nf = mpl_tab_num_flds(dca);
   for (k = 1; k <= nf; k++)
   {     switch (mpl_tab_get_type(dca, k))
      {  case 'N':
            len += 20;
            break;
         case 'S':
            len += db_escaped_string_length(mpl_tab_get_str(dca, k));
            len += 2;
            break;
              default:
                        xassert(dca != dca);
         }
   }
   query = xmalloc( (len + 1 ) * sizeof(char) );
   query[0] = 0x00;
   for (k = 1, part = strtok (template, "?"); (part != NULL);
      part = strtok (NULL, "?"), k++)
   {
      if (k > nf) break;
      strcat( query, part );
      switch (mpl_tab_get_type(dca, k))
      {  case 'N':
#if 0 /* 02/XI-2010 by xypron */
            sprintf(num, "%-18g",mpl_tab_get_num(dca, k));
#else
            sprintf(num, "%.*g", DBL_DIG, mpl_tab_get_num(dca, k));
#endif
            strcat( query, num );
            break;
         case 'S':
            strcat( query, "'");
            db_escape_string( query + strlen(query),
               mpl_tab_get_str(dca, k) );
            strcat( query, "'");
            break;
              default:
                        xassert(dca != dca);
         }
   }
   if (part != NULL)
      strcat(query, part);
   if (dl_mysql_query(sql->con, query))
   {
      xprintf("db_mysql_write: Query\n\"%s\"\nfailed.\n", query);
      xprintf("%s\n",dl_mysql_error(sql->con));
      xfree(query);
      xfree(template);
      return 1;
      }

   xfree(query);
   xfree(template);
   return 0;
   }

int db_mysql_close(TABDCA *dca, void *link)
{
   struct db_mysql *sql;

   sql = (struct db_mysql *) link;
   xassert(sql != NULL);
   dl_mysql_close(sql->con);
   if ( sql->mode == 'W' )
      xfree(sql->query);
   xfree(sql);
   dca->link = NULL;
   return 0;
}

#endif

/* eof */
