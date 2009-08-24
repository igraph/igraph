/* 
   IGraph library Java interface.
   Copyright (C) 2007  Tamas Nepusz <ntamas@rmki.kfki.hu>
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

#include "net_sf_igraph_VertexSet.h"

/************************** STATIC VARIABLES ***************************/

static jclass net_sf_igraph_VertexSet_class;
static jmethodID net_sf_igraph_VertexSet_getIdArray_mid;
static jmethodID net_sf_igraph_VertexSet_getTypeHint_mid;

/************************ CONVERSION ROUTINES **************************/

/*
 * Initializes the locally cached field IDs
 */
jint Java_net_sf_igraph_VertexSet_OnLoad(JNIEnv *env) {
  jclass cls;

  cls = (*env)->FindClass(env, JAVA_PACKAGE_PREFIX "/VertexSet");
  if (cls == 0) return JNI_ERR;
  net_sf_igraph_VertexSet_class = (*env)->NewWeakGlobalRef(env, cls);
  if (net_sf_igraph_VertexSet_class == 0) return JNI_ERR;

  net_sf_igraph_VertexSet_getIdArray_mid =
	  (*env)->GetMethodID(env, cls, "getIdArray", "()[J");
  if (net_sf_igraph_VertexSet_getIdArray_mid == 0) return JNI_ERR;

  net_sf_igraph_VertexSet_getTypeHint_mid =
	  (*env)->GetMethodID(env, cls, "getTypeHint", "()I");
  if (net_sf_igraph_VertexSet_getTypeHint_mid == 0) return JNI_ERR;

  return JNI_OK;
}

/*
 * Releases the weak references held
 */
void Java_net_sf_igraph_VertexSet_OnUnload(JNIEnv *env) {
  /*
  (*env)->DeleteWeakGlobalRef(env, net_sf_igraph_VertexSet_class);
  net_sf_igraph_VertexSet_class = 0;
  */
}

/*
 * Converts a Java VertexSet to an igraph_vs_t
 * @return:  zero if everything went fine, 1 if a null pointer was passed
 */
int Java_net_sf_igraph_VertexSet_to_igraph_vs(JNIEnv *env, jobject jobj, igraph_vs_t *result) {
  jint typeHint;
  jobject idArray;

  if (jobj == 0) {
    IGRAPH_CHECK(igraph_vs_all(result));
	return IGRAPH_SUCCESS;
  }

  typeHint = (*env)->CallIntMethod(env, jobj, net_sf_igraph_VertexSet_getTypeHint_mid);
  if (typeHint != 1 && typeHint != 2) {
    IGRAPH_CHECK(igraph_vs_all(result));
    return IGRAPH_SUCCESS;
  }
  
  idArray = (*env)->CallObjectMethod(env, jobj, net_sf_igraph_VertexSet_getIdArray_mid);
  if ((*env)->ExceptionCheck(env)) {
	return IGRAPH_EINVAL;
  }

  if (typeHint == 1) {
    /* Single vertex */
	jlong id[1];
	(*env)->GetLongArrayRegion(env, idArray, 0, 1, id);
	IGRAPH_CHECK(igraph_vs_1(result, (igraph_integer_t)id[0]));
  } else if (typeHint == 2) {
    /* List of vertices */
	jlong* ids;
	igraph_vector_t vec;
	long i, n;

	ids = (*env)->GetLongArrayElements(env, idArray, 0);
	n = (*env)->GetArrayLength(env, idArray);

	IGRAPH_VECTOR_INIT_FINALLY(&vec, n);
	for (i = 0; i < n; i++)
		VECTOR(vec)[i] = ids[i];
	IGRAPH_CHECK(igraph_vs_vector_copy(result, &vec));
	igraph_vector_destroy(&vec);
	IGRAPH_FINALLY_CLEAN(1);

	(*env)->ReleaseLongArrayElements(env, idArray, ids, JNI_ABORT);
  }

  (*env)->DeleteLocalRef(env, idArray);

  return IGRAPH_SUCCESS;
}


