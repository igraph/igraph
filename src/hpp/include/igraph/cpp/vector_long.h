/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_VECTOR_LONG_H
#define IGRAPHPP_VECTOR_LONG_H

#include <cstring>
#include <igraph/cpp/error.h>
#include <igraph/cpp/types.h>

namespace igraph {

/// C++-style wrapper around an igraph_vector_long_t
class VectorLong {
public:
    typedef long int& reference;
    typedef long int& const_reference;
    typedef long int* iterator;
    typedef const long int* const_iterator;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef long int value_type;
    typedef long int* pointer;
    typedef const long int* const_pointer;

private:
    /// The igraph_vector_long_t instance encapsulated by the wrapper
    /**
     * When the wrapper owns its contents, m_pVector points to this vector, otherwise
     * m_pVector points elsewhere and this vector is not used.
     */
    igraph_vector_long_t m_vector;

    /// Pointer to the igraph_vector_long_t instance encapsulated by the wrapper.
    /**
     * When the wrapper owns its contents, this pointer points to m_vector,
     * otherwise it points to somewhere else.
     */
    igraph_vector_long_t *m_pVector;

    /// Whether we own the given igraph_vector_long_t instance or not
    bool m_owner;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs a vector
    explicit VectorLong(size_type length = 0, pointer data = 0)
        : m_pVector(&m_vector), m_owner(true) {
        if (data) {
            IGRAPH_TRY(igraph_vector_long_init_copy(m_pVector, data, length));
        } else
            IGRAPH_TRY(igraph_vector_long_init(m_pVector, length));
    }

    /// Constructs a wrapper that takes ownership of the given vector
    VectorLong(igraph_vector_long_t vector) : m_vector(vector), m_pVector(&m_vector), m_owner(true) {}

    /// Constructs a wrapper that refers to the same vector as the given igraph_vector_long_t* pointer
    VectorLong(igraph_vector_long_t* vector, bool own=false) : m_vector(*vector), m_owner(own) {
        m_pVector = own ? &m_vector : vector;
    }

    /// Copy constructor
    VectorLong(const VectorLong& other) : m_pVector(&m_vector), m_owner(true) {
        IGRAPH_TRY(igraph_vector_long_copy(m_pVector, other.m_pVector));
    }

    /// Constructor from STL container
    template <typename InputIterator>
    VectorLong(InputIterator first, InputIterator last) {
        IGRAPH_TRY(igraph_vector_long_init(m_pVector, 0));
        while (first != last) {
            this->push_back(*first);
            ++first;
        }
    }

    /// Destroys the vector
    ~VectorLong() {
        if (m_owner)
            igraph_vector_long_destroy(m_pVector);
    }

    /******************/
    /* Static methods */
    /******************/

    /// Constructs a vector containing a sequence
    VectorLong Seq(value_type from, value_type to) {
        igraph_vector_long_t vec;
        IGRAPH_TRY(igraph_vector_long_init_seq(&vec, from, to));
        return VectorLong(vec);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Appends another vector to this one
    void append(const VectorLong& from) {
        IGRAPH_TRY(igraph_vector_long_append(m_pVector, from.m_pVector));
    }

    /// Returns an iterator pointing to the first element of the vector
    iterator begin() {
        return &(VECTOR(*m_pVector)[0]);
    }

    /// Returns an iterator pointing to the first element of the vector (const)
    const_iterator begin() const {
        return &(VECTOR(*m_pVector)[0]);
    }

    /// Returns the last element of the vector
    reference back() {
        return VECTOR(*m_pVector)[size()-1];
    }

    /// Returns the last element of the vector (const variant)
    const_reference back() const {
        return VECTOR(*m_pVector)[size()-1];
    }

    /// Finds an element using binary search in a sorted vector
    /**
     * Returns \c true if the element is found, \c false otherwise.
     * If \p pos is not null, it will be assigned to the position
     * where the element was found, or where it should be inserted
     * if the element was not found to keep the vector sorted.
     */
    bool binsearch(value_type what, long int *pos = 0) {
        return igraph_vector_long_binsearch(m_pVector, what, pos);
    }

    /// Removes all the elements from the vector
    /**
     * Note that this function sets the length of the vector to zero, but it
     * does not allocate previously used memory -- it is still kept in case
     * the vector has to grow again.
     */
    void clear() {
        igraph_vector_long_clear(m_pVector);
    }

    /// Returns whether a given element is in the vector, using linear search
    bool contains(value_type e) {
        return igraph_vector_long_contains(m_pVector, e);
    }

    /// Returns a pointer to the internal igraph_vector_long_t
    igraph_vector_long_t* c_vector() {
        return m_pVector;
    }

    /// Returns a const pointer to the internal igraph_vector_long_t
    const igraph_vector_long_t* c_vector() const {
        return m_pVector;
    }

    /// Returns whether the vector is empty
    bool empty() const {
        return igraph_vector_long_empty(m_pVector);
    }

    /// Returns an iterator pointing after the last element of the vector
    iterator end() {
        return &(VECTOR(*m_pVector)[size()]);
    }

    /// Returns an iterator pointing after the last element of the vector (const)
    const_iterator end() const {
        return &(VECTOR(*m_pVector)[size()]);
    }

    /// Fills the vector with the given item
    void fill(value_type element) {
        igraph_vector_long_fill(m_pVector, element);
    }

    /// Returns the first element of the vector
    reference front() {
        return VECTOR(*m_pVector)[0];
    }

    /// Returns the first element of the vector (const variant)
    const_reference front() const {
        return VECTOR(*m_pVector)[0];
    }

    /// Inserts an element into the vector at a given index
    void insert(long int index, value_type e) {
        IGRAPH_TRY(igraph_vector_long_insert(m_pVector, index, e));
    }

    /// Returns the minimum element of the vector
    value_type min() const {
        return igraph_vector_long_min(m_pVector);
    }

    /// Returns the maximum element of the vector
    value_type max() const {
        return igraph_vector_long_max(m_pVector);
    }

    /// Returns the maximum absolute difference between two vectors
    value_type maxdifference(const VectorLong& other) const {
        return igraph_vector_long_maxdifference(m_pVector, other.m_pVector);
    }

    /// Pops an element from the end of the vector
    value_type pop_back() {
        return igraph_vector_long_pop_back(m_pVector);
    }

    /// Prints the vector to the standard output
    void print() const {
        igraph_vector_long_print(m_pVector);
    }

    /// Adds a new element to the end of the vector
    void push_back(value_type e) {
        IGRAPH_TRY(igraph_vector_long_push_back(m_pVector, e));
    }

    /// Removes an element at the given index from the vector
    void remove(long int index) {
        igraph_vector_long_remove(m_pVector, index);
    }

    /// Removes a section of the vector
    void remove_section(long int from, long int to) {
        igraph_vector_long_remove_section(m_pVector, from, to);
    }

    /// Reserves space for the given number of elements in the vector
    void reserve(size_type newsize) {
        IGRAPH_TRY(igraph_vector_long_reserve(m_pVector, newsize));
    }

    /// Resizes the vector
    void resize(size_type newsize) {
        IGRAPH_TRY(igraph_vector_long_resize(m_pVector, newsize));
    }

    /// Reverses the vector in-place
    void reverse() {
        IGRAPH_TRY(igraph_vector_long_reverse(m_pVector));
    }

    /// Searches the vector for a given element from the given position
    bool search(long int from, value_type what, long int* pos = 0) {
        return igraph_vector_long_search(m_pVector, from, what, pos);
    }

    /// Shuffles the vector in-place
    void shuffle() {
        IGRAPH_TRY(igraph_vector_long_shuffle(m_pVector));
    }

    /// Returns the size of the vector
    size_type size() const {
        return igraph_vector_long_size(m_pVector);
    }

    /// Sorts the elements of the vector into ascending order
    void sort() {
        igraph_vector_long_sort(m_pVector);
    }

    /// Returns the sum of the elements of the vector
    value_type sum() const {
        return igraph_vector_long_sum(m_pVector);
    }

    /// Swaps the elements of this vector with another one if they are of equal length
    void swap(VectorLong& v2) {
        IGRAPH_TRY(igraph_vector_long_swap(m_pVector, v2.m_pVector));
    }

    /*************/
    /* Operators */
    /*************/

    /// Assignment operator: copies the given vector to this one
    VectorLong& operator=(const VectorLong& other) {
        IGRAPH_TRY(igraph_vector_long_update(m_pVector, other.m_pVector));
        return *this;
    }

    /// Assignment operator: copies the given array to this vector
    /**
     * It is assumed that the array has the required size.
     */
    VectorLong& operator=(const pointer other) {
        memcpy(m_pVector->stor_begin, other, size());
        return *this;
    }

    /// Equality check: returns true if the two vectors are equal
    bool operator==(const VectorLong& other) const {
        return igraph_vector_long_is_equal(m_pVector, other.m_pVector);
    }

    /// Nonequality check: returns true if the two vectors are not equal
    bool operator!=(const VectorLong& other) const {
        return !((*this) == other);
    }

    /// Returns the element with the given index
    reference operator[](long int index) {
        return VECTOR(*m_pVector)[index];
    }

    /// Returns the element with the given index (const variant)
    value_type operator[](long int index) const {
        return VECTOR(*m_pVector)[index];
    }

    /// In-place addition of a constant
    VectorLong& operator+=(value_type plus) {
        igraph_vector_long_add_constant(m_pVector, plus);
        return *this;
    }

    /// In-place addition of a vector
    VectorLong& operator+=(const VectorLong& v2) {
        IGRAPH_TRY(igraph_vector_long_add(m_pVector, v2.m_pVector));
        return *this;
    }

    /// Addition of a scalar to a vector
    VectorLong operator+(value_type plus) const {
        VectorLong result(*this);
        result += plus;
        return result;
    }

    /// In-place subtraction of a constant
    VectorLong& operator-=(value_type minus) {
        igraph_vector_long_add_constant(m_pVector, -minus);
        return *this;
    }

    /// In-place subtraction of a vector
    VectorLong& operator-=(const VectorLong& v2) {
        IGRAPH_TRY(igraph_vector_long_sub(m_pVector, v2.m_pVector));
        return *this;
    }

    /// Negation of a vector
    VectorLong operator-(void) const {
        return (*this) * -1;
    }

    /// Subtraction of a scalar from a vector
    VectorLong operator-(value_type minus) const {
        VectorLong result(*this);
        result -= minus;
        return result;
    }

    /// In-place multiplication by a constant
    VectorLong& operator*=(const real_t by) {
        igraph_vector_long_scale(m_pVector, by);
        return *this;
    }

    /// In-place division by a constant
    VectorLong& operator/=(real_t by) {
        igraph_vector_long_scale(m_pVector, 1.0 / by);
        return *this;
    }

    /// Multiplication by a constant
    VectorLong operator*(real_t by) const {
        VectorLong result(*this);
        result *= by;
        return result;
    }

    /// Division by a constant
    VectorLong operator/(real_t by) const {
        VectorLong result(*this);
        result /= by;
        return result;
    }

    /// Vector-vector scalar product
    long int operator*(const VectorLong& vector) const;
};

/// Addition of a constant to a vector from the left
VectorLong operator+(integer_t plus, const VectorLong& vector);

/// Subtraction of a constant from a vector from the left
VectorLong operator-(integer_t minus, const VectorLong& vector);

/// Multiplication of a vector by a constant from the left
VectorLong operator*(real_t by, const VectorLong& vector);

}       // end of namespaces

#endif  // IGRAPHPP_VECTOR_LONG_H

