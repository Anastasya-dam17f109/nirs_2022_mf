#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/rayleigh.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <atlimage.h>
#include <atlstr.h.>

using namespace boost::math;
using  namespace std;

const double pi = boost::math::constants::pi<double>();


struct target {
	int x;
	int y;
	int size;
	int brightness;
	int mix_type;
};

enum mix_type {
	NORMAL,
	RAYLEIGH,
	LOGNORMAL
};

class mix_img_obj
{
	shared_ptr<shared_ptr<double[]>[]> mixture_image;
	shared_ptr<double[]> re_mix_shift ;
	shared_ptr<double[]> re_mix_scale ;

	shared_ptr<wstring[]> mask_list;

	unique_ptr<wstring[]> load_mask_list;
	vector<int>           load_mask_list_idx;
	unique_ptr<wstring[]> class_list;
	wstring filename_load_image;

	mix_type mixture_type;
	string filename_gen_image = "D:\\generated_image.txt";
	CString item_type = L"item0";
	shared_ptr<target[]> targs;

	std::ofstream out;
	unsigned image_len_x = 32;
	unsigned image_len_y = 32;
	unsigned class_amount = 1;
	unsigned amount_trg = 1;
	unsigned min_targ_size = 5;
	unsigned backg_size = 32;


public:
	mix_img_obj() {};
	mix_img_obj(int img_size, mix_type mix_t, int amount_targets, int classes);
	mix_img_obj(string file_name);
	
	void     img_generator();
	void	load_from_bitmap();
	void     print_results();
	int      mean(shared_ptr<shared_ptr<double[]>[]>);

	string   get_filename()     {return filename_gen_image;}
	unsigned get_min_targ_size(){return min_targ_size;}
	double   get_bask_shift()   {return re_mix_shift[0];}
	double   get_bask_scale()   {return re_mix_scale[0];}
	mix_type get_mixture_type() {return mixture_type;}
	int      get_class_amount() {return class_amount;}
	CString  get_item_type()    {return item_type;}
	shared_ptr<target[]>  get_targets() { return targs; }
	shared_ptr<wstring[]> get_mask_list() {return mask_list; }
	std::pair<int, int>   get_image_len() {return std::pair<int, int>(image_len_x, image_len_y); }
	shared_ptr<double[]>  get_shift()     {return re_mix_shift;}
	shared_ptr<double[]>  get_scale()     {return re_mix_scale;}
	shared_ptr<shared_ptr<double[]>[]> get_image() { return mixture_image; }

	~mix_img_obj() {};
};





