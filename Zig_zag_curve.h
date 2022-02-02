#pragma once
#include "Basic_curve.h"


using namespace std;
class Zig_zag_curve: public Basic_curve
{
	int curve_type = 0;
public:
	Zig_zag_curve(int n_ord, int type);
	virtual void get_points_for_curve();
	virtual void reverse_curve() { reverse(std::begin(points), std::end(points)); }
	virtual void transpose_curve() { cout << "Basic class" << endl; }
	~Zig_zag_curve() {};
};

