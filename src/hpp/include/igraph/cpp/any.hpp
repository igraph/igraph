/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ANY_HPP
#define IGRAPHPP_ANY_HPP

#include <algorithm>
#include <iostream>
#include <typeinfo>

namespace igraph {

/// Class that can be used to hold <em>any</em> type of value
/**
 * This is similar to Boost's <tt>boost::any</tt> and it tries to
 * be API-compatible.
 */
class any {

private:
	/// Abstract base class for object containers
	class container {
	public:
		virtual ~container() {}
		/// Returns the type of the object stored in this container
		virtual const std::type_info& type() const = 0;
        /// Creates another container containing exactly the same object
        virtual container* clone() const = 0;
	};

	/// Concrete class template for object containers
	template <typename ObjectType>
	class container_impl : public container {
	public:
		/// The object being held here
		ObjectType content;

    private:
		/// Unimplemented assignment operator to prevent copying
		container_impl& operator=(const container_impl&);

	public:
		/// Creates a container that stores the given object
		container_impl(const ObjectType& value) : content(value) {}

		/// Returns the type of the object stored
		virtual const std::type_info& type() const {
			return typeid(ObjectType);
		}

		/// Creates another container storing exactly the same object
		virtual container* clone() const {
			return new container_impl(content);
		}
	};

private:
	container *content;

public:
	/// Constructs an object that does not hold anything
	any() : content(0) {}

	/// Constructs an object that stores the given value
	template <typename ObjectType>
	any(const ObjectType& obj) : content(new container_impl<ObjectType>(obj))
	{}

	/// Copy constructor
	any(const any& other) : content(other.content ? other.content->clone() : 0)
	{}

	/// Destructor
	~any() { delete content; }

	/// Checks if the object is holding anything
	bool empty() const { return content != 0; }

	/// Returns the type of the encapsulated object
	const std::type_info& type() const {
		return content ? content->type() : typeid(void);
	}

	/// Assignment operator overloading
	template <typename ObjectType>
	any& operator=(const ObjectType& obj) {
		any dummy(obj);
		std::swap(dummy.content, this->content);
		return (*this);
	}

	/// Assignment operator for "any" objects
	any& operator=(any obj) {
		std::swap(obj.content, this->content);
		return (*this);
	}
	
	/// Swaps the contents of two "any" objects
	any& swap(any& other) {
		std::swap(this->content, other.content);
		return (*this);
	}

	/// Retrieves a pointer to the object stored
	template <typename ObjectType>
	ObjectType* as() const {
		if (type() != typeid(ObjectType)) throw std::bad_cast();
		return &static_cast<container_impl<ObjectType>* >(content)->content;
	}

	/// Retrieves a reference to the object stored
	template <typename ObjectType>
	ObjectType& as() {
		if (type() != typeid(ObjectType)) throw std::bad_cast();
		return static_cast<container_impl<ObjectType>* >(content)->content;
	}
};

}           // end of namespaces

#endif      // IGRAPHPP_ANY_HPP
