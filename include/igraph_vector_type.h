/*
   igraph library.
   Copyright (C) 2013-2025  The igraph development team <igraph@igraph.org>

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

/**
 * Vector, dealing with arrays efficiently.
 * \ingroup types
 */

typedef struct TYPE(igraph_vector) {
    BASE* stor_begin;
    BASE* stor_end;
    BASE* end;
} TYPE(igraph_vector);
