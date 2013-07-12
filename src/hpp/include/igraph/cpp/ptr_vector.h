/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_PTR_VECTOR_H
#define IGRAPHPP_PTR_VECTOR_H

#include <igraph/cpp/error.h>
#include <igraph/cpp/types.h>

namespace igraph {

/// C++-style wrapper around an igraph_vector_ptr_t
template <typename T=void*>
class PtrVector {
public:
    typedef void*& reference;
    typedef T* const const_iterator;
    typedef T* iterator;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef T* pointer;

private:
    /// The igraph_vector_ptr_t instance encapsulated by the wrapper
    igraph_vector_ptr_t m_vector_ptr;

    /// Whether we own the given igraph_vector_ptr_t instance or not
    bool m_owner;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs a pointer vector
    explicit PtrVector(long int length = 0) : m_owner(true) {
        IGRAPH_TRY(igraph_vector_ptr_init(&m_vector_ptr, length));
    }

    /// Constructs a wrapper that wraps the given igraph_vector_ptr_t instance
    /**
     * \param  own  specifies whether we own the given \c igraph_vector_ptr_t instance
     *              or not. If own is false, the encapsulated \c igraph_vector_ptr_t
     *              will not be destroyed when the \c PtrVector wrapper is destroyed.
     */
    PtrVector(igraph_vector_ptr_t vector_ptr, bool own=true)
        : m_vector_ptr(vector_ptr), m_owner(own) {}

    /// Constructs a wrapper that refers to the same vector as the given igraph_vector_ptr_t* pointer
    PtrVector(igraph_vector_ptr_t* vector, bool own=false) : m_vector_ptr(*vector), m_owner(own) {}

    /// Copy constructor
    PtrVector(const PtrVector& other) : m_owner(true) {
        IGRAPH_TRY(igraph_vector_ptr_copy(&m_vector_ptr, &other.m_vector_ptr));
    }

    /// Destroys the vector
    ~PtrVector() {
        if (m_owner)
            igraph_vector_ptr_destroy(&m_vector_ptr);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Returns an iterator pointing to the first element of the vector
    iterator begin() {
        return reinterpret_cast<iterator>(&VECTOR(m_vector_ptr)[0]);
    }

    /// Returns an iterator pointing to the first element of the vector
    const_iterator begin() const {
        return reinterpret_cast<iterator>(&VECTOR(m_vector_ptr)[0]);
    }

    /// Returns the last element of the pointer vector
    reference back() {
        return VECTOR(m_vector_ptr)[size()-1];
    }

    /// Removes all the elements from the vector
    /**
     * Note that this function sets the length of the vector to zero, but it
     * does not allocate previously used memory -- it is still kept in case
     * the vector has to grow again.
     */
    void clear() {
        igraph_vector_ptr_clear(&m_vector_ptr);
    }

    /// Returns a pointer to the internal igraph_vector_ptr_t
    igraph_vector_ptr_t* c_vector_ptr() {
        return &m_vector_ptr;
    }

    /// Returns a const pointer to the internal igraph_vector_ptr_t
    const igraph_vector_ptr_t* c_vector_ptr() const {
        return &m_vector_ptr;
    }

    /// Returns whether the vector is empty
    bool empty() const {
        return size() == 0;
    }

    /// Returns an iterator pointing after the last element of the vector
    iterator end() {
        return begin() + size();
    }

    /// Returns an iterator pointing after the last element of the vector (const)
    const_iterator end() const {
        return begin() + size();
    }

    /// Returns the first element of the vector
    reference front() {
        return VECTOR(m_vector_ptr)[0];
    }

    /// Returns an element of the vector cast to the given pointer type
    value_type get(long int index) {
        return static_cast<T>(VECTOR(m_vector_ptr)[index]);
    }

    /// Returns an element of the vector cast to the given pointer type (const)
    value_type const get(int index) const {
        return static_cast<T const>(VECTOR(m_vector_ptr)[index]);
    }

    /// Adds a new element to the end of the vector
    void push_back(void* e) {
        IGRAPH_TRY(igraph_vector_ptr_push_back(&m_vector_ptr, e));
    }

    /// Resizes the vector
    void resize(long int newsize) {
        IGRAPH_TRY(igraph_vector_ptr_resize(&m_vector_ptr, newsize));
    }

    /// Returns the size of the vector
    size_t size() const {
        return igraph_vector_ptr_size(&m_vector_ptr);
    }

    /*************/
    /* Operators */
    /*************/

    /// Assignment operator; intentionally left unimplemented.
    PtrVector& operator=(const PtrVector& other);

    /// Returns the element with the given index
    reference operator[](long int index) {
        return VECTOR(m_vector_ptr)[index];
    }
};

}       // end of namespaces

#endif  // IGRAPHPP_VECTOR_H

