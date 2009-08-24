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

#include "jni_utils.h"
#include "config.h"
#include <igraph/igraph.h>
#include <string.h>          /* strlen */

/************************** STATIC VARIABLES ***************************/

static JavaVM *jvm;

/*********************** INITIALIZER FUNCTION **************************/

JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *this_jvm, void* reserved) {
  JNIEnv *env;

  jvm = this_jvm;
  if ((*jvm)->GetEnv(jvm, (void**)&env, JNI_VERSION_1_2)) return JNI_ERR;

  if (Java_net_sf_igraph_Graph_OnLoad(env) == JNI_ERR) return JNI_ERR;

  if (Java_net_sf_igraph_Connectedness_OnLoad(env) == JNI_ERR) return JNI_ERR;
  if (Java_net_sf_igraph_NeighborMode_OnLoad(env) == JNI_ERR) return JNI_ERR;
  if (Java_net_sf_igraph_StarMode_OnLoad(env) == JNI_ERR) return JNI_ERR;

  if (Java_net_sf_igraph_VertexSet_OnLoad(env) == JNI_ERR) return JNI_ERR;

  return JNI_VERSION_1_2;
}

JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *this_jvm, void* reserved) {
  JNIEnv *env;

  if ((*jvm)->GetEnv(jvm, (void**)&env, JNI_VERSION_1_2)) return;

  Java_net_sf_igraph_Graph_OnUnload(env);

  Java_net_sf_igraph_Connectedness_OnUnload(env);
  Java_net_sf_igraph_NeighborMode_OnUnload(env);
  Java_net_sf_igraph_StarMode_OnUnload(env);

  Java_net_sf_igraph_VertexSet_OnUnload(env);
}

/************************ AUXILIARY FUNCTIONS **************************/

/*
 * Returns the environment of the current thread using the cached JVM
 */
JNIEnv *JNU_GetEnv() {
  JNIEnv *env;
  (*jvm)->GetEnv(jvm, (void**)&env, JNI_VERSION_1_2);
  return env;
}

/*
 * Throws an exception by exception class name
 * Adapted from http://java.sun.com/docs/books/jni/html/exceptions.html#26050
 */
void JNU_ThrowByName(JNIEnv *env, const char *name, const char *msg) {
  jclass cls = (*env)->FindClass(env, name);
  if (cls != 0)
	(*env)->ThrowNew(env, cls, msg);
  (*env)->DeleteLocalRef(env, cls);
}

/********** THINGS TO DO BEFORE ENTERING & AFTER LEAVING C LAYER ***********/

static igraph_error_handler_t *Java_igraph_old_error_handler;
static igraph_warning_handler_t *Java_igraph_old_warning_handler;

void Java_igraph_error_handler(const char *reason, const char *file,
  int line, int igraph_errno) {
  JNIEnv *env = JNU_GetEnv();
  char msg[8192], *p;
  IGRAPH_FINALLY_FREE();

  if ((*env)->ExceptionCheck(env)) {
	/* We already have an exception, keep that and return */
	return;
  }

  if (strlen(reason) > 2 && reason[0] == '_' && reason[1] == '_') {
	/* Special case: throwing a Java exception by name */
	JNU_ThrowByName(env, reason+2, "");
  } else {
    snprintf(msg, 8192, "%s, %s at %s:%i", reason, igraph_strerror(igraph_errno),
      file, line);
    JNU_ThrowByName(env, JAVA_PACKAGE_PREFIX "/CoreException", msg);
  }
}

void Java_igraph_before() {
  Java_igraph_old_error_handler=igraph_set_error_handler(Java_igraph_error_handler);
  Java_igraph_old_warning_handler=igraph_set_warning_handler(0);
}

void Java_igraph_after() {
  igraph_set_error_handler(Java_igraph_old_error_handler);
  igraph_set_warning_handler(Java_igraph_old_warning_handler);
}

