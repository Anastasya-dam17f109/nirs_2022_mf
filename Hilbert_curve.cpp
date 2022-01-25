#include "pch.h"
#include "Hilbert_curve.h"

// конструктор

Hilbert_curve::Hilbert_curve(int n_ord) :
	Basic_curve(n_ord){}

// создание точки кривой заданного порядка

Point Hilbert_curve::from_d(int step) {
	Point p = { 0, 0 };
	bool rx, ry;
	int t = step;
	for (int s = 1; s < n_order; s <<= 1) {
		rx = ((t & 2) != 0);
		ry = (((t ^ (rx ? 1 : 0)) & 1) != 0);
		p.rot(s, rx, ry);
		p.x += (rx ? s : 0);
		p.y += (ry ? s : 0);
		t >>= 2;
	}
	return p;
}

// получение массива точек кривой

void Hilbert_curve::get_points_for_curve() {
	for (int d = 0; d < n_order * n_order; ++d) 
		points.push_back(from_d(d));
}

// отрисовка точек  Гильбертовой кривой

std::vector<std::string> Hilbert_curve::draw_curve() {
	int row, col, delta_X, delta_Y;
	auto canvas = new char *[n_order];
	for (size_t i = 0; i < n_order; i++) {
		canvas[i] = new char[n_order * 3 - 2];
		memset(canvas[i], ' ', n_order * 3 - 2);
	}

	for (int i = 1; i < points.size(); i++) {
		auto last_point = points[i - 1];
		auto cur_point = points[i];
		delta_X = cur_point.x - last_point.x;
		delta_Y = cur_point.y - last_point.y;
		if (delta_X == 0) {
			// vertical line
			row = max(cur_point.y, last_point.y);
			col = cur_point.x * 3;
			canvas[row][col] = '|';
		}
		else {
			// horizontal line
			row = cur_point.y;
			col = min(cur_point.x, last_point.x) * 3 + 1;
			canvas[row][col] = '_';
			canvas[row][col + 1] = '_';
		}
	}

	vector<std::string> lines;
	for (size_t i = 0; i < n_order; i++) {
		string temp;
		temp.assign(canvas[i], n_order * 3 - 2);
		lines.push_back(temp);
	}

	for (size_t i = 0; i < n_order; i++)
		delete[] canvas[i];
	delete[] canvas;
	return lines;
}

// перегрузка оператора присваивания

Hilbert_curve Hilbert_curve::operator= (const Hilbert_curve drob){
	n_order = drob.n_order;
	points  = drob.points;
	return *this;
}


