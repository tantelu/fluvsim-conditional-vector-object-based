#pragma once
#include "basedefine.h"
#include "spline.h"

using namespace tk;

class __declspec(dllexport) curveinterp
{
private:
	spline* line = NULL;

public:
	curveinterp(const vector<Coordinate>& ctrls);
	~curveinterp();

	Coordinate get_dot(double t);
};

