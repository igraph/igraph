/* platform.c
 *
 * Copyright (C) 2014 Tamas Nepusz
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

#include "platform.h"

#ifdef _MSC_VER

inline double _plfit_fmin(double a, double b) {
	return (a < b) ? a : b;
}

inline double _plfit_round(double x) {
	return floor(x+0.5);
}

#endif

/* Dummy function to prevent a warning when compiling with Clang - the file
 * would contain no symbols */
void _plfit_i_unused() {}
