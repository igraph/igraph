#include "unit_limiter.h"

double unit_limiter(double vUnitDouble)
{
	double result = vUnitDouble;
	if (result < 0.0)
		result = 0.0;
	else if (result > 1.0)
		result = 1.0;
	return result;
}
