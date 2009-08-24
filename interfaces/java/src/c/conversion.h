/* 
   IGraph library Java interface.
   Copyright (C) 2007-2009  Tamas Nepusz <tamas@cs.rhul.ac.uk>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA 

*/

/*

ATTENTION: This is a highly experimental, proof-of-concept Java interface.
Its main purpose was to convince me that it can be done in finite time :)
The interface is highly incomplete, at the time of writing even some
essential functions (e.g. addEdges) are missing. Since I don't use Java
intensively, chances are that this interface gets finished only if there
is substantial demand for it and/or someone takes the time to send patches
or finish it completely.

*/

#ifndef _Included_net_sf_igraph_conversion
#define _Included_net_sf_igraph_conversion

#ifdef __cplusplus
extern "C" {
#endif

#include <jni.h>
#include <igraph/igraph.h>

/************************ CONVERSION ROUTINES **************************/

/* Conversion between jobject and igraph_t */
int Java_jobject_to_igraph(JNIEnv *env, jobject jobj, igraph_t** gptr);
jobject Java_igraph_to_new_jobject(JNIEnv *env, igraph_t* gptr, jclass cls);

/* Conversion between jdoubleArray and igraph_vector_t */
int Java_jdoubleArray_to_igraph_vector(JNIEnv *env, jdoubleArray array, igraph_vector_t* vector);
jdoubleArray Java_igraph_vector_to_new_jdoubleArray(JNIEnv *env, igraph_vector_t* vector);

#ifdef __cplusplus
}
#endif
#endif

