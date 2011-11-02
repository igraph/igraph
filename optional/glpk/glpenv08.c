/* glpenv08.c (shared library support) */

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glpenv.h"

/* GNU version ********************************************************/

#if defined(HAVE_LTDL)

#include <ltdl.h>

void *xdlopen(const char *module)
{     void *h = NULL;
      if (lt_dlinit() != 0)
      {  lib_err_msg(lt_dlerror());
         goto done;
      }
      h = lt_dlopen(module);
      if (h == NULL)
      {  lib_err_msg(lt_dlerror());
         if (lt_dlexit() != 0)
            xerror("xdlopen: %s\n", lt_dlerror());
      }
done: return h;
}

void *xdlsym(void *h, const char *symbol)
{     void *ptr;
      xassert(h != NULL);
      ptr = lt_dlsym(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: %s\n", symbol, lt_dlerror());
      return ptr;
}

void xdlclose(void *h)
{     xassert(h != NULL);
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
{     void *h;
      h = dlopen(module, RTLD_NOW);
      if (h == NULL)
         lib_err_msg(dlerror());
      return h;
}

void *xdlsym(void *h, const char *symbol)
{     void *ptr;
      xassert(h != NULL);
      ptr = dlsym(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: %s\n", symbol, dlerror());
      return ptr;
}

void xdlclose(void *h)
{     xassert(h != NULL);
      if (dlclose(h) != 0)
         xerror("xdlclose: %s\n", dlerror());
      return;
}

/* Windows version ****************************************************/

#elif defined(__WOE__)

#include <windows.h>

void *xdlopen(const char *module)
{     void *h;
      h = LoadLibrary(module);
      if (h == NULL)
      {  char msg[20];
         sprintf(msg, "Error %d", GetLastError());
         lib_err_msg(msg);
      }
      return h;
}

void *xdlsym(void *h, const char *symbol)
{     void *ptr;
      xassert(h != NULL);
      ptr = GetProcAddress(h, symbol);
      if (ptr == NULL)
         xerror("xdlsym: %s: Error %d\n", symbol, GetLastError());
      return ptr;
}

void xdlclose(void *h)
{     xassert(h != NULL);
      if (!FreeLibrary(h))
         xerror("xdlclose: Error %d\n", GetLastError());
      return;
}

/* NULL version *******************************************************/

#else

void *xdlopen(const char *module)
{     xassert(module == module);
      lib_err_msg("Shared libraries not supported");
      return NULL;
}

void *xdlsym(void *h, const char *symbol)
{     xassert(h != h);
      xassert(symbol != symbol);
      return NULL;
}

void xdlclose(void *h)
{     xassert(h != h);
      return;
}

#endif

/* eof */
