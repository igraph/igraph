/** Vector.h
 */

#ifndef VECTOR_H
#define VECTOR_H

#include "Point.h"

namespace igraph {

class Vector
{
public:
	Vector();
	Vector(const Point& vStartPoint, const Point& vEndPoint);
	Vector(double vI, double vJ, double vK);
	~Vector();

	Vector Normalize() const; // returns a unit vector of this vector
	void NormalizeThis();
	void ReverseDirection();
	bool IsSameDirection(const Vector& rVector) const;
	
	void I(double vI);
	double I() const;
	void J(double vJ);
	double J() const;
	void K(double vK);
	double K() const;

	double Dot(const Vector& rVector) const; // returns the dot product of this and rVector
	Vector Cross(const Vector& rVector) const; // returns the cross product of this and rVector
	Vector operator+ (Vector vRhs) const; // returns the sum of two vectors
	Vector operator- (Vector vRhs) const; // returns the difference of two vectors
	Point operator+ (Point vRhs) const; // returns the sum of a vector and a Point 
	Vector operator* (double vRhs) const; // returns multiplication of a scalar with a vector
	Point ToPoint() const;	// converts a vector to a point

	double Magnitude() const;

private:


	double mI, mJ, mK;
};

} // namespace igraph

#endif
