/* tls.c (thread local storage) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2001-2017 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "env.h"

#ifndef TLS
static void *tls = NULL;
#else
static TLS void *tls = NULL;
/* this option allows running multiple independent instances of GLPK in
 * different threads of a multi-threaded application, in which case the
 * variable tls should be placed in the Thread Local Storage (TLS);
 * it is assumed that the macro TLS is previously defined to something
 * like '__thread', '_Thread_local', etc. */
#endif

/***********************************************************************
*  NAME
*
*  tls_set_ptr - store global pointer in TLS
*
*  SYNOPSIS
*
*  #include "env.h"
*  void tls_set_ptr(void *ptr);
*
*  DESCRIPTION
*
*  The routine tls_set_ptr stores a pointer specified by the parameter
*  ptr in the Thread Local Storage (TLS). */

void tls_set_ptr(void *ptr)
{     tls = ptr;
      return;
}

/***********************************************************************
*  NAME
*
*  tls_get_ptr - retrieve global pointer from TLS
*
*  SYNOPSIS
*
*  #include "env.h"
*  void *tls_get_ptr(void);
*
*  RETURNS
*
*  The routine tls_get_ptr returns a pointer previously stored by the
*  routine tls_set_ptr. If the latter has not been called yet, NULL is
*  returned. */

void *tls_get_ptr(void)
{     void *ptr;
      ptr = tls;
      return ptr;
}

/**********************************************************************/

#ifdef __WOE__

/*** Author: Heinrich Schuchardt <xypron.glpk@gmx.de> ***/

#pragma comment(lib, "user32.lib")

#include <windows.h>

#define VISTA 0x06

/* This is the main entry point of the DLL. */

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID
      lpvReserved)
{     DWORD version;
      DWORD major_version;
#ifdef TLS
      switch (fdwReason)
      {  case DLL_PROCESS_ATTACH:
         /* @TODO:
          * GetVersion is deprecated but the version help functions are
          * not available in Visual Studio 2010. So lets use it until
          * we remove the outdated Build files. */
         version = GetVersion();
         major_version = version & 0xff;
         if (major_version < VISTA)
         {  MessageBoxA(NULL,
               "The GLPK library called by this application is configur"
               "ed to use thread local storage which is not fully suppo"
               "rted by your version of Microsoft Windows.\n\n"
               "Microsoft Windows Vista or a later version of Windows i"
               "s required to run this application.",
               "GLPK", MB_OK | MB_ICONERROR);
            return FALSE;
         }
         break;
      }
#endif /* TLS */
      return TRUE;
}

#endif /* __WOE__ */

/* eof */
