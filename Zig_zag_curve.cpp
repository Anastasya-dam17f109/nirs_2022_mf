#include "pch.h"
#include "Zig_zag_curve.h"


Zig_zag_curve::Zig_zag_curve(int n_ord, int type): Basic_curve(n_ord){
	curve_type = type;
}

//создание масива точек кривой-зигзага в зависимости от типа зигзага (тип обуславливает  граничные и начальные условия)

void Zig_zag_curve::get_points_for_curve() {

	Point cur_point = {0, 0 };
	int f_step, s_step, f_rearward, s_rearward, f_front, s_front, pre_front;
	int * f_coord = &cur_point.x; int * s_coord = &cur_point.x;
	switch (curve_type) {
		case 0: {
			f_step = 1;
			s_step = 1;
			f_rearward = 0;
			s_rearward = 0;
			f_front = n_order - 1;
			s_front = n_order - 1;
			pre_front = n_order - 2;
			f_coord = &cur_point.x;
			s_coord = &cur_point.y;
		}
		break;
		case 1:{
			f_step = 1;
			s_step = 1;
			f_rearward = 0;
			s_rearward = 0;
			f_front = n_order - 1;
			s_front = n_order - 1;
			pre_front = n_order - 2;
			f_coord = &cur_point.y;
			s_coord = &cur_point.x;
		}
		break;

		case 2: {
			cur_point.x = n_order - 1;
			cur_point.y = 0;
			f_step = -1;
			s_step = 1;
			f_rearward = n_order - 1;
			s_rearward = 0;
			f_front = 0;
			s_front = n_order - 1;
			
			pre_front = 1;
			f_coord = &cur_point.x;
			s_coord = &cur_point.y;
		}
		break;
		case 3: {
			cur_point.x = n_order - 1;
			cur_point.y = 0;
			f_step = 1;
			s_step = -1;
			f_rearward = 0;
			s_rearward = n_order - 1;
			f_front = n_order - 1;
			s_front = 0;
			pre_front = n_order - 1-1;
			f_coord = &cur_point.y;
			s_coord = &cur_point.x;
		}
		break;
	}
	points.push_back(cur_point);
	if (n_order != 1) {
		while (*f_coord != pre_front) {
			(*f_coord)+= f_step;
			points.push_back(cur_point);
			while ((*f_coord) != f_rearward) {
				(*f_coord) -= f_step;
				(*s_coord) += s_step;
				points.push_back(cur_point);
			}
			(*s_coord)+= s_step;
			points.push_back(cur_point);
			while (*s_coord != s_rearward) {
				(*s_coord)-= s_step;
				(*f_coord)+= f_step;
				points.push_back(cur_point);
			}

		}
		(*f_coord)+= f_step;
		(*s_coord) = s_rearward -s_step;

		while (((*f_coord) != pre_front) && ((*s_coord) != s_front)) {
			(*s_coord)+= s_step;
			points.push_back(cur_point);
			while ((*s_coord) != s_front) {
				(*f_coord)-= f_step;
				(*s_coord)+= s_step;
				points.push_back(cur_point);
			}
			(*f_coord)+= f_step;
			points.push_back(cur_point);
			while ((*f_coord) != f_front) {
				(*s_coord)-= s_step;
				(*f_coord)+= f_step;
				points.push_back(cur_point);
			}
		}
	}
}

