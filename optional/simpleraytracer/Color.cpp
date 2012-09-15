#include "Color.h"
#include "unit_limiter.h"

Color::Color()
{
}

Color::Color(double vRed, double vGreen, double vBlue)
{
	Red(vRed);
	Green(vGreen);
	Blue(vBlue);
}

Color::~Color()
{
}

// returns multiplication of a scalar with this vector
Color Color::operator* (double vRhs) const
{
	return Color(mRed*vRhs, mGreen*vRhs, mBlue*vRhs);
}

// returns the addition of this color with another color
Color Color::operator+ (const Color& vRhs) const
{
	return Color(Red()+vRhs.Red(),Green()+vRhs.Green(),Blue()+vRhs.Blue());
}

void Color::Red(double vRed)
{
	mRed = unit_limiter(vRed);
}
double Color::Red() const
{
	return mRed;
}
void Color::Green(double vGreen)
{
	mGreen = unit_limiter(vGreen);

}
double Color::Green() const
{
	return mGreen;
}
void Color::Blue(double vBlue)
{
	mBlue = unit_limiter(vBlue);
}
double Color::Blue() const
{
	return mBlue;
}

unsigned char Color::RedByte() const
{
	return ByteValue(mRed);
}
unsigned char Color::GreenByte() const
{
	return ByteValue(mGreen);
}
unsigned char Color::BlueByte() const
{
	return ByteValue(mBlue);
}

unsigned char Color::ByteValue(double vZeroToOne) const
{
	return (unsigned char)(vZeroToOne*255.0);
}
