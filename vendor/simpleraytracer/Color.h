/** Color.h
 */

#ifndef COLOR_H
#define COLOR_H

namespace igraph {

class Color
{
public:
	Color();
	Color(double vRed, double vGreen, double vBlue, 
	      double vTransparent=1.0);
	~Color();
	
	Color operator* (double vRhs) const; // returns multiplication of a scalar with a vector
	Color operator+ (const Color& vRhs) const; // returns the addition of this color with another color

	void Red(double vRed);
	double Red() const;
	void Green(double vGreen);
	double Green() const;
	void Blue(double vBlue);
	double Blue() const;
	void Transparent(double vTransparent);
	double Transparent() const;

	unsigned char RedByte() const;
	unsigned char GreenByte() const;
	unsigned char BlueByte() const;
	unsigned char TransparentByte() const;
private:
	unsigned char ByteValue(double vZeroToOne) const;
	double mRed, mGreen, mBlue, mTransparent;
};

} // namespace igraph

#endif
