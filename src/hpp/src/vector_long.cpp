/* vim:set ts=4 sw=4 sts=4 et: */

#include <cassert>
#include <igraph/cpp/vector_long.h>

namespace igraph {

long int VectorLong::operator*(const VectorLong& vector) const {
    VectorLong::value_type result = VectorLong::value_type();

    assert(vector.size() == size());

    VectorLong::const_iterator it1 = begin(), it2 = vector.begin(), it3 = end();
    while (it1 != it3) {
        result += (*it1) * (*it2);
        it1++; it2++;
    }
    return result;
}

VectorLong operator+(integer_t plus, const VectorLong& vector) {
    return vector * plus;
}

VectorLong operator-(integer_t minus, const VectorLong& vector) {
    return (-vector) + minus;
}

VectorLong operator*(real_t by, const VectorLong& vector) {
    return vector * by;
}

}         // end of namespaces
