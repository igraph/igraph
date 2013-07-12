/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_MATRIX_H
#define IGRAPHPP_MATRIX_H

#include <cstring>
#include <igraph/igraph_blas.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

/// C++-style wrapper around an igraph_matrix_t
class Matrix {
public:
    typedef igraph_real_t* iterator;
    typedef const igraph_real_t* const_iterator;

private:
    /// The igraph_matrix_t instance encapsulated by the wrapper
    igraph_matrix_t m_matrix;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an empty matrix
    explicit Matrix(long nrow = 0, long ncol = 0) {
        IGRAPH_TRY(igraph_matrix_init(&m_matrix, nrow, ncol));
    }

    /// Copy constructor
    Matrix(const Matrix& other) {
        IGRAPH_TRY(igraph_matrix_copy(&m_matrix, &other.m_matrix));
    }

    /// Destroys the matrix
    ~Matrix() {
        igraph_matrix_destroy(&m_matrix);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Returns an iterator pointing to the first element of the matrix
    iterator begin() {
        return &(MATRIX(m_matrix, 0, 0));
    }

    /// Returns an iterator pointing to the first element of the matrix (const)
    const_iterator begin() const {
        return &(MATRIX(m_matrix, 0, 0));
    }

    /// Returns a const pointer to the internal igraph_matrix_t
    const igraph_matrix_t* c_matrix() const {
        return &m_matrix;
    }

    /// Returns a pointer to the internal igraph_matrix_t
    igraph_matrix_t* c_matrix() {
        return &m_matrix;
    }

    /// Returns the sums of the columns
    Vector colsum() const {
        Vector result;
        colsum(result);
        return result;
    }

    /// Returns the sums of the columns
    void colsum(Vector& result) const {
        IGRAPH_TRY(igraph_matrix_colsum(&m_matrix, result.c_vector()));
    }

    /// Returns an iterator pointing after the last element of the matrix
    iterator end() {
        return &(MATRIX(m_matrix, 0, 0)) + size();
    }

    /// Returns an iterator pointing after the last element of the matrix (const)
    const_iterator end() const {
        return &(MATRIX(m_matrix, 0, 0)) + size();
    }

    /// Fills the matrix with the given element
    void fill(igraph_real_t element) {
        igraph_matrix_fill(&m_matrix, element);
    }

    /// Returns the given column of the matrix as a vector
    Vector getCol(long int index) const {
        Vector result(nrow());
        IGRAPH_TRY(igraph_matrix_get_col(&m_matrix, result.c_vector(), index));
        return result;
    }

    /// Returns the given column of the matrix as a vector
    void getCol(long int index, Vector& result) const {
        result.resize(nrow());
        IGRAPH_TRY(igraph_matrix_get_col(&m_matrix, result.c_vector(), index));
    }

    /// Returns the given row of the matrix as a vector
    Vector getRow(long int index) const {
        Vector result(ncol());
        IGRAPH_TRY(igraph_matrix_get_row(&m_matrix, result.c_vector(), index));
        return result;
    }

    /// Returns the given row of the matrix as a vector
    void getRow(long int index, Vector& result) const {
        result.resize(ncol());
        IGRAPH_TRY(igraph_matrix_get_row(&m_matrix, result.c_vector(), index));
    }

    /// Returns whether the matrix is symmetric
    bool isSymmetric() const {
        return igraph_matrix_is_symmetric(&m_matrix);
    }

    /// Returns the minimum element of the matrix
    igraph_real_t min() const {
        return igraph_matrix_min(&m_matrix);
    }

    /// Returns the maximum element of the matrix
    igraph_real_t max() const {
        return igraph_matrix_max(&m_matrix);
    }

    /// Returns the maximum absolute difference between two matrices
    igraph_real_t maxdifference(const Matrix& other) const {
        return igraph_matrix_maxdifference(&m_matrix, &other.m_matrix);
    }

    /// Returns the number of columns of the matrix
    size_t ncol() const {
        return igraph_matrix_ncol(&m_matrix);
    }

    /// Returns the number of rows of the matrix
    size_t nrow() const {
        return igraph_matrix_nrow(&m_matrix);
    }

    /// Prints the matrix to the standard output
    void print() const {
        igraph_matrix_print(&m_matrix);
    }

    /// Resizes the matrix
    void resize(long nrow, long ncol) {
        IGRAPH_TRY(igraph_matrix_resize(&m_matrix, nrow, ncol));
    }

    /// Returns the sums of the rows
    Vector rowsum() const {
        Vector result;
        rowsum(result);
        return result;
    }

    /// Returns the sums of the rows
    void rowsum(Vector& result) const {
        IGRAPH_TRY(igraph_matrix_rowsum(&m_matrix, result.c_vector()));
    }

    /// Returns the size of the matrix
    size_t size() const {
        return igraph_matrix_size(&m_matrix);
    }

    /// Returns the sum of the elements of the matrix
    igraph_real_t sum() const {
        return igraph_matrix_sum(&m_matrix);
    }

    /*************/
    /* Operators */
    /*************/

    /// Assignment operator: copies the given matrix to this one
    Matrix& operator=(const Matrix& other) {
        IGRAPH_TRY(igraph_matrix_update(&m_matrix, &other.m_matrix));
        return *this;
    }

    /// Assignment operator: copies the given array
    /**
     * It is assumed that the array has the required size and it is
     * stored in the same ordering as the matrix.
     */
    Matrix& operator=(const igraph_real_t* other) {
        memcpy(m_matrix.data.stor_begin, other, sizeof(igraph_real_t) * size());
        return *this;
    }

    /// Equality check: returns true if the two matrices are equal
    bool operator==(const Matrix& other) const {
        return igraph_matrix_is_equal(&m_matrix, &other.m_matrix);
    }

    /// Inequality check: returns true if the two matrices are not equal
    bool operator!=(const Matrix& other) const {
        return !((*this) == other);
    }

    /// Returns the element with the given row and column indices
    igraph_real_t& operator()(long int ri, long int ci) {
        return MATRIX(m_matrix, ri, ci);
    }

    /// Returns the element with the given row and column indices (const variant)
    igraph_real_t& operator()(long int ri, long int ci) const {
        return MATRIX(m_matrix, ri, ci);
    }

    /// Multiplication by vector from right
    Vector operator*(const Vector& v) const {
        Vector result(v.size());
        igraph_blas_dgemv(0, 1, &m_matrix, v.c_vector(),
                0, result.c_vector());
        return result;
    }

    /// In-place addition with constant
    Matrix& operator+=(igraph_real_t plus) {
        igraph_matrix_add_constant(&m_matrix, plus);
        return *this;
    }

    /// In-place addition with matrix
    Matrix& operator+=(const Matrix& other) {
        IGRAPH_TRY(igraph_matrix_add(&m_matrix, &other.m_matrix));
        return *this;
    }

    /// In-place subtraction
    Matrix& operator-=(const Matrix& other) {
        IGRAPH_TRY(igraph_matrix_sub(&m_matrix, &other.m_matrix));
        return *this;
    }

    /// In-place scaling
    Matrix& operator*=(igraph_real_t by) {
        igraph_matrix_scale(&m_matrix, by);
        return *this;
    }

    /// In-place division
    Matrix& operator/=(igraph_real_t by) {
        igraph_matrix_scale(&m_matrix, 1.0/by);
        return *this;
    }
};

}       // end of namespaces

#endif  // IGRAPHPP_VECTOR_H

