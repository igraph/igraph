/*
   igraph library.
   Copyright (C) 2016-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#undef IGRAPH_BEGIN_C_DECLS
#undef IGRAPH_END_C_DECLS
#ifdef __cplusplus
    #define IGRAPH_BEGIN_C_DECLS extern "C" {
    #define IGRAPH_END_C_DECLS }
#else
    #define IGRAPH_BEGIN_C_DECLS /* empty */
    #define IGRAPH_END_C_DECLS /* empty */
#endif

/* This is to eliminate gcc warnings about unused parameters */
#define IGRAPH_UNUSED(x) (void)(x)

/* The pure function attribute of GCC-compatible compilers indicates
 * that the function does not have side-effects, i.e. it does not
 * modify global memory. This enables additional compiler optimizations
 * such as common subexpression elimination.
 *
 * The const attribute is similar but with much more stringent requirements.
 * The function must also not read global memory. Generally, const functions
 * should not take pointers, and must compute the return value solely based
 * on their input.
 */
#ifdef __GNUC__
#define IGRAPH_FUNCATTR_PURE __attribute__((__pure__))
#define IGRAPH_FUNCATTR_CONST __attribute__((__const__))
#else
#define IGRAPH_FUNCATTR_PURE
#define IGRAPH_FUNCATTR_CONST
#endif

/* IGRAPH_ASSUME() provides hints to the compiler about conditions
 * that are true yet the compiler cannot deduce. Use with great care.
 * Assuming a condition that is not actually true leads to undefined behaviour. */
#if defined(__clang__)
   /* For Clang, see https://clang.llvm.org/docs/LanguageExtensions.html */
#  if __has_builtin(__builtin_assume)
#    define IGRAPH_ASSUME(expr) __builtin_assume(expr)
#  else
#    define IGRAPH_ASSUME(expr) /* empty */
#  endif
#elif defined(__GNUC__) && !defined(__ICC)
   /* Introduced in GCC 4.5, https://gcc.gnu.org/gcc-4.5/changes.html */
#  define IGRAPH_ASSUME(expr) do { if (expr) {} else { __builtin_unreachable(); } } while (0)
#elif defined(_MSC_VER) || defined(__ICC)
#  define IGRAPH_ASSUME(expr) __assume(expr)
#else
#  define IGRAPH_ASSUME(expr) /* empty */
#endif

/* IGRAPH_I_STRINGIFY(X) evaluates X and converts the result to a string. */
#define IGRAPH_I_STRINGIFY_I(X) #X
#define IGRAPH_I_STRINGIFY(X) IGRAPH_I_STRINGIFY_I(X)

/* Include the definition of macros controlling symbol visibility */
#include "igraph_export.h"

/* Used instead of IGRAPH_EXPORT with functions that need to be tested,
 * but are not part of the public API. */
#define IGRAPH_PRIVATE_EXPORT IGRAPH_EXPORT

/* Marks experimental functions. If IGRAPH_WARN_EXPERIMENTAL is defined to a
 * nonzero value, supported compilers will emit a warning when these functions
 * are used. */
#if defined(__GNUC__) && IGRAPH_WARN_EXPERIMENTAL
#define IGRAPH_EXPERIMENTAL __attribute__((__warning__("Experimental function.")))
#else
#define IGRAPH_EXPERIMENTAL /* empty */
#endif
