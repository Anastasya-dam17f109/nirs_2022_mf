#ifndef COMBINE_CURVE_H
#define COMBINE_CURVE_H
#include "Basic_curve.h"
#include "Zig_zag_curve.h"
#include "Hilbert_curve.h"

using namespace std;
class combine_curve: public Basic_curve
{
public:
	combine_curve(int n_ord, int type, bool  flag);
	virtual void get_points_for_curve()
	{ 
		switch (curve_type)
		{
		case 0:
		{
			Zig_zag_curve c1(n_order, 0, false);
			c1.get_points_for_curve();
			for (auto& i : c1.get_points())
				points.push_back(i);
			Zig_zag_curve c2(n_order, 0, true);
			c2.get_points_for_curve();
			for (auto& i : c2.get_points())
				points.push_back(i);

			Hilbert_curve c3(n_order, 0, false);
			c3.get_points_for_curve();
			for (auto& i : c3.get_points())
				points.push_back(i);

			Zig_zag_curve c4(n_order, 2, false);
			c4.get_points_for_curve();
			for (auto& i : c4.get_points())
				points.push_back(i);

			Zig_zag_curve c5(n_order, 2, true);
			c5.get_points_for_curve();
			for (auto& i : c5.get_points())
				points.push_back(i);

			Hilbert_curve c6(n_order, 0, true);
			c6.get_points_for_curve();
			for (auto& i : c6.get_points())
				points.push_back(i);
		}
		break;
		case 1:
		{
			Hilbert_curve c6(n_order, 0, true);
			c6.get_points_for_curve();
			for (auto& i : c6.get_points())
				points.push_back(i);
			Zig_zag_curve c1(n_order, 0, false);
			c1.get_points_for_curve();
			for (auto& i : c1.get_points())
				points.push_back(i);
			Zig_zag_curve c2(n_order, 0, true);
			c2.get_points_for_curve();
			for (auto& i : c2.get_points())
				points.push_back(i);
			Hilbert_curve c3(n_order, 0, false);
			c3.get_points_for_curve();
			for (auto& i : c3.get_points())
				points.push_back(i);

			Zig_zag_curve c4(n_order, 2, false);
			c4.get_points_for_curve();
			for (auto& i : c4.get_points())
				points.push_back(i);

			Zig_zag_curve c5(n_order, 2, true);
			c5.get_points_for_curve();
			for (auto& i : c5.get_points())
				points.push_back(i);

		}
		break;
		}
		
	}
	
	~combine_curve() {};
};


#endif
