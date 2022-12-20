//#include "pch.h"
#include "combine_curve.h"
//#include "Zig_zag_curve.h"
//#include "Hilbert_curve.h"

combine_curve::combine_curve(int n_ord, int type, bool flag): Basic_curve(n_ord, type, flag){
}

//создание масива точек кривой-зигзага в зависимости от типа зигзага (тип обуславливает  граничные и начальные условия)

//void combine_curve::get_points_for_curve() {

	/*Zig_zag_curve c1(n_order, 0, false));
	c1.get_points_for_curve();
	points.push_back(c1.get_points());
	Zig_zag_curve c2(n_order, 0, true));
	c2.get_points_for_curve();
	points.push_back(c2.get_points());

	Hilbert_curve c3(n_order, 0, false);
	c3.get_points_for_curve();
	points.push_back(c3.get_points());

	Zig_zag_curve c4(n_order, 2, false));
	c4.get_points_for_curve();
	points.push_back(c4.get_points());

	Zig_zag_curve c5(n_order, 2, true));
	c5.get_points_for_curve();
	points.push_back(c5.get_points());

	Hilbert_curve c6(n_order, 0, true);
	c6.get_points_for_curve();
	points.push_back(c6.get_points());*/
//	int a = 5;
//}

