#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

/* In igraph 0.8, we use DECLDIR only with MSVC, not other compilers on Windows. */
#undef DECLDIR
#if defined (_MSC_VER)
    #ifdef IGRAPH_EXPORTS
        #define DECLDIR __declspec(dllexport)
    #elif defined(IGRAPH_STATIC)
        #define DECLDIR /**/
    #else
        #define DECLDIR __declspec(dllimport)
    #endif
#else
    #define DECLDIR /**/
#endif
