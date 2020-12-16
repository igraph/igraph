#include "prpack_result.h"
#include <cstdlib>
using namespace prpack;

prpack_result::prpack_result() {
    x = NULL;
}

prpack_result::~prpack_result() {
    delete[] x;
}

