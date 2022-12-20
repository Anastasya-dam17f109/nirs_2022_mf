#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

// точка кривой Гильберта

struct Point {
	int x, y;

	//rotate/flip a quadrant appropriately
	void rot(int n, bool rx, bool ry) {
		if (!ry) {
			if (rx) {
				x = (n - 1) - x;
				y = (n - 1) - y;
			}
			swap(x, y);
		}
	}
};

class Basic_curve{
public:
	int n_order = 1;
	int curve_type = 0;
	bool reverse_flag = false;
	vector<Point> points;
	Basic_curve() {}
	Basic_curve(int n_ord, int type, bool flag):n_order(n_ord), curve_type(type), reverse_flag(flag){};
	virtual void get_points_for_curve() { cout << "Basic class" << endl; }
	virtual void reverse_curve() { cout << "Basic class" << endl; }
	virtual void transpose_curve() { cout << "Basic class" << endl; }
	virtual void reflect_curve(bool x_flag) { cout << "Basic class" << endl; }
	vector<Point>& get_points() { return  points; }
	~Basic_curve() {};
};

