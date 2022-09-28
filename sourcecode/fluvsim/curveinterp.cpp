#include "curveinterp.h"

curveinterp::curveinterp(const vector<Coordinate>& ctrls)
{
	if (ctrls.size() < 3) {
		throw exception("The number of points must be greater than 3");
	}
	vector<double>x;
	vector<double>y;
	for (size_t i = 0; i < ctrls.size(); i++)
	{
		x.push_back(ctrls[i].x);
		y.push_back(ctrls[i].y);
	}
	line = new spline(x, y);
}

curveinterp::~curveinterp()
{
	if (line != NULL) {
		delete line;
		line = NULL;
	}
}

Coordinate curveinterp::get_dot(double x)
{
	auto y = (*line)(x);
	return Coordinate(x, y);
}
