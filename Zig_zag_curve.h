#pragma once
#include "Basic_curve.h"


using namespace std;
class Zig_zag_curve: public Basic_curve
{
	int curve_type = 0;
public:
	Zig_zag_curve(int n_ord, int type);
	void get_points_for_curve();
	~Zig_zag_curve() {};
};

