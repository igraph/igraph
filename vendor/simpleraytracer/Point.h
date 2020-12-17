/**
	this is a simple generic class representing a 3d point with a name.
	it also defines the PointList type, which is a linked list of Points
*/

#ifndef POINT_H
#define POINT_H

#include <list>
using namespace std;

namespace igraph {

class Point
{
public:
	Point(); // creates a point at the origin with name 0
	Point(double vX, double vY, double vZ, int vName);
	Point(double vX, double vY, double vZ);
	~Point();

	double X() const;
	void X(double vX);
	double Y() const;
	void Y(double vY);
	double Z() const;
	void Z(double vZ);

	int Name() const;
	void Name(int vName);
	double Distance(const Point& rPoint) const;

	bool operator==(const Point& vRhs) const;

private:
	double mX, mY, mZ;
	int mName;	
};

typedef list<Point> PointList;
typedef list<Point>::iterator PointListIterator;

} // namespace igraph

#endif
