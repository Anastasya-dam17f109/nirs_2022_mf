#ifndef HELIX_CURVE_H
#define HELIX_CURVE_H
#include "Basic_curve.h"

class Helix_curve: public Basic_curve
{
public:
    Helix_curve(){}
    Helix_curve(int n_ord, int type, bool flag):Basic_curve( n_ord, type, flag){}
    virtual void get_points_for_curve();
    virtual void reverse_curve(){reverse(std::begin(points), std::end(points));}


};

#endif // HELIX_CURVE_H
