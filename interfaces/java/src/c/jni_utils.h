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

#ifndef _Included_net_sf_igraph_jni_utils
#define _Included_net_sf_igraph_jni_utils

#ifdef __cplusplus
extern "C" {
#endif

#include <jni.h>
#include "config.h"

/*********************** INITIALIZER FUNCTION **************************/

JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *this_jvm, void* reserved);
JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *this_jvm, void* reserved);

/************************ AUXILIARY FUNCTIONS **************************/

/// Returns the environment of the current thread using the cached JVM
JNIEnv *JNU_GetEnv();

/// Throws an exception by exception class name
void JNU_ThrowByName(JNIEnv *env, const char *name, const char *msg);

/********** THINGS TO DO BEFORE ENTERING & AFTER LEAVING C LAYER ***********/

void Java_igraph_error_handler(const char *reason, const char *file, int line, int igraph_errno);
void Java_igraph_before();
void Java_igraph_after();

#ifdef __cplusplus
}
#endif
#endif

