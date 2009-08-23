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

#include "conversion.h"

extern jfieldID net_sf_igraph_Graph_handle_fid;
extern jmethodID net_sf_igraph_Graph_constructor_mid;

/***** Conversion between jobject and igraph_t */

/*
 * Converts a Java jobject to an igraph_t* if appropriate by casting
 * its private handle field to an igraph_t*
 *
 * @return: 0 if everything was OK, 1 otherwise
 */
int Java_jobject_to_igraph(JNIEnv *env, jobject jobj, igraph_t** gptr) {
  *gptr = (igraph_t*)(uintptr_t)((*env)->GetLongField(env, jobj, net_sf_igraph_Graph_handle_fid));
  return (*gptr == 0 ? 1 : 0);
}

/*
 * Converts an igraph_t* to a new Java Graph object
 * @return: the new Java Graph object or NULL if there was an error
 */
jobject Java_igraph_to_new_jobject(JNIEnv *env, igraph_t* gptr, jclass cls) {
  /* Construct the object */
  jobject result;
  result = (*env)->NewObject(env, cls, net_sf_igraph_Graph_constructor_mid, gptr);
  return result;
}

/***** Conversion between jdoubleArray and igraph_vector_t */

/**
 * Converts a Java double[] to an igraph_vector_t* object.
 * The igraph_vector_t* that's passed in must be uninitialized.
 * @return: zero if everything was OK, an igraph error code otherwise
 */
int Java_jdoubleArray_to_igraph_vector(JNIEnv *env, jdoubleArray array, igraph_vector_t* vector) {
	jsize i, n = (*env)->GetArrayLength(env, array);
	jdouble* elements = (*env)->GetDoubleArrayElements(env, array, 0);

	IGRAPH_CHECK(igraph_vector_init(vector, n));
	for (i=0; i < n; i++)
		VECTOR(*vector)[i] = elements[i];

	(*env)->ReleaseDoubleArrayElements(env, array, elements, JNI_ABORT);

	return 0;
}

/**
 * Converts an igraph_vector_t* to a new Java double[] object
 * @return: the new Java double array or NULL if there was an error
 */
jdoubleArray Java_igraph_vector_to_new_jdoubleArray(JNIEnv *env, igraph_vector_t* vector) {
	long i, n;
	jdoubleArray result;
	jdouble* elements;

	n = igraph_vector_size(vector);
	result = (*env)->NewDoubleArray(env, n);
	elements = (*env)->GetDoubleArrayElements(env, result, 0);

	for (i=0; i < n; i++) {
		elements[i] = VECTOR(*vector)[i];
	}

	(*env)->ReleaseDoubleArrayElements(env, result, elements, 0);

	return result;
}

