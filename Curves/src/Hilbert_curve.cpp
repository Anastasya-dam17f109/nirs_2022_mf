#include "pch.h"
#include "Hilbert_curve.h"

// конструктор

Hilbert_curve::Hilbert_curve(int n_ord, int type, bool flag) :
	Basic_curve(n_ord, type, flag){
	
}

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
	switch (curve_type) {
		case 0:{}
		break;
		case 1: {
			// поворот на 90против часовой - как транспонирование
			transpose_curve();
		}
		break;
		case 2: {
			// отражение относительно гризонтальной прямой
			reflect_curve(false);
		}
		break;
		case 3: {
			// поворот на 90по часовой - как транспонирование
			transpose_curve();
			reflect_curve(true);
		}
		break;
	}
	if (reverse_flag)
		reverse_curve();
}

// поворот-транспонирование

void Hilbert_curve::transpose_curve() {
	int buf;
	for (auto &point : points) {
		buf = point.x;
		point.x = point.y;
		point.y = buf;
	}
}


// зеркальное отражение относительно оси(горизонталь пока)

void Hilbert_curve::reflect_curve(bool x_flag) {
	int buf_dist;
	int mid = n_order / 2;
	if (x_flag) {
		for (auto &point : points) 
			point.x = n_order - 1 - point.x;
	}
	else {
		for (auto &point : points)
			point.y = n_order - 1 - point.y;
	}
	/*std::vector<std::string> lines = draw_curve();
		for (auto &line : lines) {
			std::cout << line.c_str() << '\n';
		}*/
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


