/* dlsup.c (dynamic linking support) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2008-2013 Free Software Foundation, Inc.
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

/* GNU version ********************************************************/

#if defined(HAVE_LTDL)

#include <ltdl.h>

void *xdlopen(const char *module)
{     /* open dynamically linked library */
      void *h = NULL;
      if (lt_dlinit() != 0)
      {  put_err_msg(lt_dlerror());
         goto done;
      }
      h = lt_dlopen(module);
      if (h == NULL)
      {  put_err_msg(lt_dlerror());
         if (lt_dlexit() != 0)
            xerror("xdlopen: %s\n", lt_dlerror());
      }
done: return h;
}

void *xdlsym(void *h, const char *symbol)
{     /* obtain address of symbol from dynamically linked library */
      void *ptr;
      xassert(h != NULL);
      ptr = lt_dlsym(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: %s\n", symbol, lt_dlerror());
      return ptr;
}

void xdlclose(void *h)
{     /* close dynamically linked library */
      xassert(h != NULL);
      if (lt_dlclose(h) != 0)
         xerror("xdlclose: %s\n", lt_dlerror());
      if (lt_dlexit() != 0)
         xerror("xdlclose: %s\n", lt_dlerror());
      return;
}

/* POSIX version ******************************************************/

#elif defined(HAVE_DLFCN)

#include <dlfcn.h>

void *xdlopen(const char *module)
{     /* open dynamically linked library */
      void *h;
      h = dlopen(module, RTLD_NOW);
      if (h == NULL)
         put_err_msg(dlerror());
      return h;
}

void *xdlsym(void *h, const char *symbol)
{     /* obtain address of symbol from dynamically linked library */
      void *ptr;
      xassert(h != NULL);
      ptr = dlsym(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: %s\n", symbol, dlerror());
      return ptr;
}

void xdlclose(void *h)
{     /* close dynamically linked library */
      xassert(h != NULL);
      if (dlclose(h) != 0)
         xerror("xdlclose: %s\n", dlerror());
      return;
}

/* MS Windows version *************************************************/

#elif defined(__WOE__)

#include <windows.h>

void *xdlopen(const char *module)
{     /* open dynamically linked library */
      void *h;
      h = LoadLibrary(module);
      if (h == NULL)
      {  char msg[20];
         sprintf(msg, "Error %d", GetLastError());
         put_err_msg(msg);
      }
      return h;
}

void *xdlsym(void *h, const char *symbol)
{     /* obtain address of symbol from dynamically linked library */
      void *ptr;
      xassert(h != NULL);
      ptr = GetProcAddress(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: Error %d\n", symbol, GetLastError());
      return ptr;
}

void xdlclose(void *h)
{     /* close dynamically linked library */
      xassert(h != NULL);
      if (!FreeLibrary(h))
         xerror("xdlclose: Error %d\n", GetLastError());
      return;
}

/* NULL version *******************************************************/

#else

void *xdlopen(const char *module)
{     /* open dynamically linked library */
      xassert(module == module);
      put_err_msg("Shared libraries not supported");
      return NULL;
}

void *xdlsym(void *h, const char *symbol)
{     /* obtain address of symbol from dynamically linked library */
      xassert(h != h);
      xassert(symbol != symbol);
      return NULL;
}

void xdlclose(void *h)
{     /* close dynamically linked library */
      xassert(h != h);
      return;
}

#endif

/* eof */
