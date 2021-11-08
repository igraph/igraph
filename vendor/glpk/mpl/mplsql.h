/* mplsql.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2003-2016 Free Software Foundation, Inc.
*  Written by Heinrich Schuchardt <xypron.glpk@gmx.de>.
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

#ifndef MPLSQL_H
#define MPLSQL_H

#define db_iodbc_open _glp_db_iodbc_open
void *db_iodbc_open(TABDCA *dca, int mode);
/* open iODBC database connection */

#define db_iodbc_read _glp_db_iodbc_read
int db_iodbc_read(TABDCA *dca, void *link);
/* read data from iODBC */

#define db_iodbc_write _glp_db_iodbc_write
int db_iodbc_write(TABDCA *dca, void *link);
/* write data to iODBC */

#define db_iodbc_close _glp_db_iodbc_close
int db_iodbc_close(TABDCA *dca, void *link);
/* close iODBC database connection */

#define db_mysql_open _glp_db_mysql_open
void *db_mysql_open(TABDCA *dca, int mode);
/* open MySQL database connection */

#define db_mysql_read _glp_db_mysql_read
int db_mysql_read(TABDCA *dca, void *link);
/* read data from MySQL */

#define db_mysql_write _glp_db_mysql_write
int db_mysql_write(TABDCA *dca, void *link);
/* write data to MySQL */

#define db_mysql_close _glp_db_mysql_close
int db_mysql_close(TABDCA *dca, void *link);
/* close MySQL database connection */

#endif

/* eof */
