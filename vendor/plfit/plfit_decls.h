/* plfit_decls.h
 *
 * Copyright (C) 2024 Tamas Nepusz
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

#ifndef PLFIT_DECLS_H
#define PLFIT_DECLS_H

#undef PLFIT_BEGIN_C_DECLS
#undef PLFIT_END_C_DECLS
#ifdef __cplusplus
    #define PLFIT_BEGIN_C_DECLS extern "C" {
    #define PLFIT_END_C_DECLS }
#else
    #define PLFIT_BEGIN_C_DECLS /* empty */
    #define PLFIT_END_C_DECLS /* empty */
#endif

#define PLFIT_EXPORT /* empty */

#endif /* PLFIT_DECLS_H */
