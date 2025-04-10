/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_constructors.h"
#include "igraph_operators.h"

#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"


igraph_error_t igraph_i_cartesian_product(igraph_t *res,
   const igraph_t *g1, const igraph_t *g2) {

   return IGRAPH_SUCCESS;  // dummy return                        
}

igraph_error_t igraph_product(igraph_t *res,
                        const igraph_t *g1, const igraph_t *g2,
                        igraph_product_t type) {

   switch (type) {
      case IGRAPH_CARTESIAN_PRODUCT:
         IGRAPH_CHECK(igraph_i_cartesian_product(res, g1, g2));
         return IGRAPH_SUCCESS;
      
      default:
         IGRAPH_ERROR("Unknown graph product type.", IGRAPH_EINVAL);
   }
}
