/* vim:set ts=4 sw=4 sts=4 et: */

#include <cassert>
#include <igraph/cpp/vector_bool.h>

namespace igraph {

long int VectorBool::operator*(const VectorBool& vector) const {
    long int result = 0;

    assert(vector.size() == size());

    VectorBool::const_iterator it1 = begin(), it2 = vector.begin(), it3 = end();
    while (it1 != it3) {
        if (*it1 && *it2)
            result++;
        it1++; it2++;
    }
    return result;
}

}         // end of namespaces
