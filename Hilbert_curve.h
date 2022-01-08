#pragma once
#include "Basic_curve.h"

class Hilbert_curve : public Basic_curve
{
public:
	Hilbert_curve(): Basic_curve() {}
	Hilbert_curve(int n_ord);
	Point from_d(int step);
	void get_points_for_curve();
	vector<string> draw_curve();
	//Hilbert_curve& operator= (const Hilbert_curve &drob);
	Hilbert_curve operator= (const Hilbert_curve drob);
	~Hilbert_curve() {};
};

