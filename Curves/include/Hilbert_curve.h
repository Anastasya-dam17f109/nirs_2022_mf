#pragma once
#include "Basic_curve.h"

class Hilbert_curve : public Basic_curve
{
	
public:
	Hilbert_curve(): Basic_curve() {}
	Hilbert_curve(int n_ord, int  type, bool flag);
	Point from_d(int step);
	virtual void get_points_for_curve();
	vector<string> draw_curve();
	//Hilbert_curve& operator= (const Hilbert_curve &drob);
	virtual void reverse_curve() {reverse(points.begin(), points.end());}
	virtual void transpose_curve(); 
	virtual void reflect_curve(bool x_flag);
	Hilbert_curve operator= (const Hilbert_curve drob);
	~Hilbert_curve() {};
};

