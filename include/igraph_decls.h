#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

#undef DECLDIR
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    #if defined (__MINGW32__) || defined (__CYGWIN32__)
        #define DECLDIR /**/
    #else
        #ifdef IGRAPH_EXPORTS
            #define DECLDIR __declspec(dllexport)
        #elif defined(IGRAPH_STATIC)
            #define DECLDIR /**/
        #else
            #define DECLDIR __declspec(dllimport)
        #endif
    #endif
#else
    #define DECLDIR /**/
#endif
