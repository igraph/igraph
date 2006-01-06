/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef IGRAPH_ERROR_H
#define IGRAPH_ERROR_H

/* This file contains the igraph error handling.
 * Most bits are taken literally from the GSL library (with the GSL_
 * prefix renamed to IGRAPH_), as i couldn't find a better way to do
 * them. */

/**
 * SECTION:errorhandling
 * @title: Error handling
 * @short_description: how <literal>igraph</literal> handles errors
 *
 * </para>
 * <refsect2><title>Error handlers</title>
 * <para>
 * If igraph runs into an error -- an invalid argument was supplied
 * to a function, or we've runned out of memory -- the control is
 * transferred to the error handler function.
 * </para><para>
 * The default error handler is #igraph_error_handler_abort which
 * prints and error message and aborts the program.
 * </para>
 * <para>
 * The igraph_set_error_handler() function can be used to set a new
 * error handler function of type #igraph_error_handler_t, see the
 * documentation of this type for details.
 * </para>
 * <para>
 * There are two other predefined error handler functions,
 * #igraph_error_handler_ignore and #igraph_error_handler_printignore, 
 * these deallocate the temporarily allocated memory (more about this
 * later) and return with the error code. The latter also prints an
 * error message. If you use these error handlers you need to take
 * care about possible errors yourself by checking the return value of
 * every <literal>igraph</literal> function.
 * </para><para>
 * Independently of the error handler installed, all functions in the
 * library do their best to leave their arguments
 * <emphasis>semantically</emphasis> unchanged if an error
 * happens. (Eg. vector_reserve() leaves its argument semantically
 * unchanged, even if it allocated some more memory for it, but
 * vector_resize() does not.) 
 * </para></refsect2>
 * <refsect2><title>Error codes</title>
 * <para>
 * Every <literal>igraph</literal> function which can fail return a
 * single integer error code. Some functions are very simple and
 * cannot run into any error, these may return other types, or
 * <type>void</type> as well. The error codes are defined by the
 * #igraph_i_error_type_t enumeration.
 * </para>
 * </refsect2>
 * <refsect2><title>Writing error handlers</title>
 * <para>
 * You can write and install error handlers simply by defining a
 * function of type #igraph_error_handler_type and calling
 * igraph_set_error_handler(). This feature is useful for interface
 * writers, as the <literal>igraph</literal> will have the chance to
 * signal errors the appropriate way, eg. the R interface defines an
 * error handler which calls the <function>error()</function>
 * function, as required by R, while the Python interface has an error
 * handler which raises an exception according to the Python way.
 * </para>
 * <para> 
 * If you want to write an error handler, your error handler should
 * call IGRAPH_FINALLY_FREE() to deallocate all temporarily memory to
 * prevent memory leaks.
 * </para>
 * </refsect2>
 * <refsect2><title>Error handling internals</title>
 * <para>
 * If an error happens, the functions in the library call the
 * IGRAPH_ERROR() macro with a textual description of the error and an
 * <literal>igraph</literal> error code. igraph_error() simply calls
 * the installed error handler. Other useful macro is IGRAPH_CHECK(),
 * this check the return value of its argument which is normally a
 * function call, and calls IGRAPH_ERROR() if it is not
 * %IGRAPH_SUCCESS.
 * </para>
 * </refsect2><refsect2><title>Deallocating memory</title>
 * <para>
 * If a function runs into an error (and the program is not aborted)
 * the error handler should deallocate all temporarily memory. This is
 * done by storing the address and destroy function of all temporary
 * object in a stack. The IGRAPH_FINALLY() function declares an object as
 * temporary by placing its address in the stack. If a function returns
 * with success it calls IGRAPH_FINALLY_CLEAN() with the
 * number of objects to remove from the stack. If an error happens
 * however, the error handler should call IGRAPH_FINALLY_FREE() to
 * deallocate each object added to the stack. This means that the
 * temporary objects allocated in the calling function (and etc.) will
 * be freed as well.
 * </para></refsect2>
 * <refsect2><title>Writing <literal>igraph</literal> functions with
 * proper error handling</title>
 * <para>
 * There are some simple rules to keep in order to have functions
 * behaving well in errorenous situations. First, check the arguments
 * of the functions and call IGRAPH_ERROR() if they are invalid. Second,
 * call IGRAPH_FINALLY() on each dinamically allocated object and call
 * IGRAPH_FINALLY_CLEAN() with the proper argument. Third, use
 * IGRAPH_CHECK on all function calls which can generate errors.
 * </para>
 * <para>
 * The size of the stack used for this bookkeeping is fixed, and
 * small. If you want to allocate several objects, write a destroy
 * function which can deallocate all of these. See the
 * <filename>adjlist.c</filename> file in the
 * <literal>igraph</literal> source for an example.
 * </para>
 * <para> 
 * For some functions these mechanisms are simply not flexible
 * enough. These functions should define their own error handlers and
 * restore the error handler before they return.
 * </para>
 * </refsect2>
 * <refsect2><title>Error handling and threads</title>
 * <para>
 * It is likely that the <literal>igraph</literal> error handling
 * method is <emphasis>not</emphasis> thread-safe, mainly beacuse of
 * the static global stack which is used to store the address of the
 * temporarily allocated objects. This issue might be addressed in a
 * later version of <literal>igraph</literal>.
 * </para>
 * </refsect2><para>
 */

/**
 * igraph_error_handler_t:
 * @reason: Textual description of the error.
 * @file: The source file in which the error is noticed.
 * @line: The number of the line in the source file which triggered
 *   the error
 * @igraph_errno: The <literal>igraph</literal> error code.
 * 
 * The type of error handler functions.
 */

typedef void igraph_error_handler_t (const char * reason, const char * file,
				     int line, int igraph_errno);

/**
 * igraph_error_handler_abort: 
 *
 * The default error handler, prints an error message and aborts the
 * program. 
 */
/**
 * igraph_error_handler_ignore: 
 *
 * This error handler frees the temporarily allocated memory and returns
 * with the error code. 
 */
/**
 * igraph_error_handler_printignore: 
 * 
 * Frees temporarily allocated memory and returns with the error
 * code. 
 */

extern igraph_error_handler_t igraph_error_handler_abort;
extern igraph_error_handler_t igraph_error_handler_ignore;
extern igraph_error_handler_t igraph_error_handler_printignore;

/**
 * igraph_set_error_handler:
 * @new_handler: The error handler function to install.
 *
 * Installs a new error handler. If called with 0, it installs the
 * default error handler (which is cirrently
 * #igraph_error_handler_abort). 
 * 
 * Returns: the old error handler function. This should be saved and
 * restored if %new_handler is not needed any more.
 */

igraph_error_handler_t *
igraph_set_error_handler(igraph_error_handler_t new_handler);

/**
 * igraph_i_error_type_t:
 * @IGRAPH_SUCCESS: The function successfully completed its task.
 * @IGRAPH_FAILURE: Something went wrong. You'll almost never meet this
 * error as normally more specific error codes are used.
 * @IGRAPH_ENOMEM: There wasn't enough memory to allocate memory on the heap.
 * @IGRAPH_PARSEERROR: A parse error was found in a file.
 * @IGRAPH_EINVAL: A parameter's value is invalid. Eg. negative number
 * was specified as the number of vertices.
 * @IGRAPH_EXISTS: A graph/vertex/edge attribute is already installed
 * with the given name.
 * @IGRAPH_EINVEVECTOR: Invalid vector of vertex ids. A vertex id is
 * either negative or bigger than the number of vertices minus one.
 * @IGRAPH_EINVVID: Invalid vertex id, negative or too big.
 * @IGRAPH_NONSQUARE: A non-square matrix was received while a square
 * matrix was expected.
 * @IGRAPH_EINVMODE: Invalid mode parameter.
 * @IGRAPH_EFILE: A file operation failed. Eg. a file doesn't exist,
 * or the user ha no rights to open it.
 * @IGRAPH_EUNFOLDINF: Attempted to unfold an infinite iterator.
 */

typedef enum {
  IGRAPH_SUCCESS    = 0,
  IGRAPH_FAILURE    = 1,
  IGRAPH_ENOMEM     = 2,
  IGRAPH_PARSEERROR = 3,
  IGRAPH_EINVAL     = 4,
  IGRAPH_EXISTS     = 5,
  IGRAPH_EINVEVECTOR= 6,
  IGRAPH_EINVVID    = 7,
  IGRAPH_NONSQUARE  = 8,
  IGRAPH_EINVMODE   = 9,
  IGRAPH_EFILE      = 10,
  IGRAPH_EUNFOLDINF = 11
} igraph_i_error_type_t;

/**
 * IGRAPH_ERROR:
 * @reason: Textual description of the error. This should be something
 * more explaning than the text associated with the error code. Eg. if
 * the error code is %IGRAPH_EINVAL, its asssociated text (see
 * igraph_strerror()) is "Invalid value" and %return should explain
 * which parameter was invalid and maybe why.
 * @igraph_errno: The <literal>igraph</literal> error code.
 * 
 * This macro is called if an error is noticed. It calles
 * igraph_error() with the proper parameters and if that returns 
 * the macro returns the "calling" function as well, with the error
 * code. If for some (suspicious) reason you want to call the error
 * handler without returning from the current function, call
 * igraph_error() directly.
 */

#define IGRAPH_ERROR(reason, igraph_errno) \
       do { \
       igraph_error (reason, __FILE__, __LINE__, igraph_errno) ; \
       return igraph_errno ; \
       } while (0)

/**
 * igraph_error:
 * @reason: Textual description of the error.
 * @file: The source file in which the error was noticed.
 * @line: The number of line in the source file which triggered the
 * error.
 * 
 * This is the function which is called (usually via the
 * IGRAPH_ERROR() macro if an error is noticed. See the discussion of
 * this macro as well.
 * 
 * Returns: the error code (if it returns)
 */

int igraph_error(const char *reason, const char *file, int line,
		 int igraph_errno);

/**
 * igraph_strerror:
 * @igraph_errno: The <literal>igraph</literal> error code.
 * 
 * This is a simple utility function, it gives a short general textual
 * description for an <literal>igraph</literal> error code.
 * 
 * Returns: pointer to the textual description of the error code.
 */

const char* igraph_strerror(const int igraph_errno);

#define IGRAPH_ERROR_SELECT_2(a,b)       ((a) != IGRAPH_SUCCESS ? (a) : ((b) != IGRAPH_SUCCESS ? (b) : IGRAPH_SUCCESS))
#define IGRAPH_ERROR_SELECT_3(a,b,c)     ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_2(b,c))
#define IGRAPH_ERROR_SELECT_4(a,b,c,d)   ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_3(b,c,d))
#define IGRAPH_ERROR_SELECT_5(a,b,c,d,e) ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_4(b,c,d,e))

/* Now comes the more conveninent error handling macro arsenal.
 * Ideas taken from exception.{h,c} by Laurent Deniau see
 * http://cern.ch/Laurent.Deniau/html/oopc/oopc.html#Exceptions for more 
 * information. We don't use the exception handling code though.  */

struct igraph_i_protectedPtr {
  int all;
  void *ptr;
  void (*func)(void*);
};

typedef void igraph_finally_func_t (void*);

void IGRAPH_FINALLY_REAL(void (*func)(void*), void* ptr);

/**
 * IGRAPH_FINALLY_CLEAN:
 * @num: The number of objects to remove from the bookkeeping stack.
 * 
 * Removes the specified number of objects from the stack of
 * temporarily allocated objects. Most often this is called just
 * before returning from a function.
 */

void IGRAPH_FINALLY_CLEAN(int num); 

/**
 * IGRAPH_FINALLY_FREE:
 *
 * Calls the destroy function for all objects in the stack of
 * temporarily allocated objects. This is usually called only from an
 * error handler. It is <emphasis>not</emphasis> appropriate to use it
 * instead of destroying each unneeded object of a function, as it
 * destroys the temporary objects of the caller function (and so on)
 * as well.
 */

void IGRAPH_FINALLY_FREE();

/**
 * IGRAPH_FINALLY:
 * @func: The address of the function which is normally called to
 * destroy the object.
 * @ptr: Pointer to the object itself.
 * 
 * This macro places the address of an object, together with the
 * address of its destructor in a stack. This stack is used if an
 * error happens to deallocate temporarily allocated objects to
 * prevent memory leaks.
 */

#define IGRAPH_FINALLY(func, ptr) \
  IGRAPH_FINALLY_REAL((igraph_finally_func_t*)(func), (ptr))

/**
 * IGRAPH_CHECK:
 * @a: An expression, usually a function call.
 * 
 * Executes the expression and checks its value. If this is not
 * %IGRAPH_SUCCESS, it calls IGRAPH_ERROR() with the value as the
 * error code. Here is an example usage:
 * <informalexample>
 * <programlisting>
 * IGRAPH_CHECK(vector_push_back(&amp;v, 100));
 * </programlisting>
 * </informalexample>
 */

#define IGRAPH_CHECK(a) do { \
                 int igraph_i_ret=(a); \
                 if (igraph_i_ret != 0) { \
                    IGRAPH_ERROR("", igraph_i_ret); \
                 } } while(0)

#endif
