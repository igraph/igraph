/* plfit_error.h
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef PLFIT_ERROR_H
#define PLFIT_ERROR_H

#include "plfit_decls.h"

PLFIT_BEGIN_C_DECLS

enum {
	PLFIT_SUCCESS  = 0,
	PLFIT_FAILURE  = 1,
	PLFIT_EINVAL   = 2,
	PLFIT_UNDRFLOW = 3,
	PLFIT_OVERFLOW = 4,
	PLFIT_ENOMEM   = 5,
	PLFIT_EMAXITER = 6
};

#if (defined(__GNUC__) && GCC_VERSION_MAJOR >= 3)
#  define PLFIT_UNLIKELY(a) __builtin_expect((a), 0)
#  define PLFIT_LIKELY(a)   __builtin_expect((a), 1)
#else
#  define PLFIT_UNLIKELY(a) a
#  define PLFIT_LIKELY(a)   a
#endif

#define PLFIT_CHECK(a) \
	do {\
		int plfit_i_ret=(a); \
		if (PLFIT_UNLIKELY(plfit_i_ret != PLFIT_SUCCESS)) {\
			return plfit_i_ret; \
		} \
	} while (0)

#define PLFIT_ERROR(reason,plfit_errno) \
	do {\
		plfit_error (reason, __FILE__, __LINE__, plfit_errno) ; \
		return plfit_errno ; \
	} while (0)

typedef void plfit_error_handler_t(const char*, const char*, int, int);

PLFIT_EXPORT extern plfit_error_handler_t plfit_error_handler_abort;
PLFIT_EXPORT extern plfit_error_handler_t plfit_error_handler_ignore;
PLFIT_EXPORT extern plfit_error_handler_t plfit_error_handler_printignore;

PLFIT_EXPORT plfit_error_handler_t* plfit_set_error_handler(plfit_error_handler_t* new_handler);

PLFIT_EXPORT void plfit_error(const char *reason, const char *file, int line, int plfit_errno);
PLFIT_EXPORT const char* plfit_strerror(const int plfit_errno);

PLFIT_END_C_DECLS

#endif /* PLFIT_ERROR_H */
