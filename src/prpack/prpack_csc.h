#ifndef PRPACK_CSC
#define PRPACK_CSC

#ifndef _MSC_VER
#  ifdef HAVE_STDINT_H
#    include <stdint.h>
#  else
#    ifdef HAVE_SYS_INT_TYPES_H
#      include <sys/int_types.h>
#    else
#      include "pstdint.h"
#    endif
#  endif
#else
typedef __int64 int64_t;
#endif

namespace prpack {

    class prpack_csc {
        public:
            int num_vs;
            int num_es;
            int* heads;
            int* tails;
    };

    class prpack_int64_csc {
        public:
            int64_t num_vs;
            int64_t num_es;
            int64_t* heads;
            int64_t* tails;
    };
};

#endif
