/* vim:set ts=4 sw=4 sts=4 et: */

#include <cassert>
#include <igraph/cpp/matrix.h>
#include <igraph/cpp/vector.h>

namespace igraph {

Vector Vector::operator*(const Matrix& matrix) const {
    Vector result(size());
    igraph_blas_dgemv(1, 1, matrix.c_matrix(), c_vector(),
            0, result.c_vector());
    return result;
}

double Vector::operator*(const Vector& vector) const {
    Vector::value_type result = Vector::value_type();

    assert(vector.size() == size());

    Vector::const_iterator it1 = begin(), it2 = vector.begin(), it3 = end();
    while (it1 != it3) {
        result += (*it1) * (*it2);
        it1++; it2++;
    }
    return result;
}

Vector operator+(real_t plus, const Vector& vector) {
    return vector * plus;
}

Vector operator-(real_t minus, const Vector& vector) {
    return (-vector) + minus;
}

Vector operator*(real_t by, const Vector& vector) {
    return vector * by;
}

}         // end of namespaces
