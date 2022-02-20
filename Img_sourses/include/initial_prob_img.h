#pragma once
#include "mix_img_obj.h"


class initial_prob_img
{
	shared_ptr<mix_img_obj> m_image;
	shared_ptr <int[]>      layer_size;
	shared_ptr <int[]>      layer_idx;
	shared_ptr <int[]>      init_layer_idx;
	int layer_amount = 1;
	double* init_prob_img;
	shared_ptr<double[]> mix_shift;
	shared_ptr<double[]> mix_scale;
	shared_ptr<double[]> cl_probs;
	unsigned image_len_x = 32;
	unsigned image_len_y = 32;
	unsigned class_amount = 1;

public:
	initial_prob_img(int img_size, mix_type mix_t, int amount_targets, int classes);
	initial_prob_img(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes);
	initial_prob_img(shared_ptr<mix_img_obj> image);
	initial_prob_img(string filename);
	void     alloc_layer_mmr();
	void     generate_init_probs_mixtures();

	int                     get_class_amount() { return class_amount; }
	std::pair<int, int>     get_image_len() { return std::pair<int, int>(image_len_x, image_len_y); }
	double*                 get_image() { return init_prob_img;}
	shared_ptr<mix_img_obj> get_m_image(){return m_image; }
	shared_ptr<double[]>    get_cl_probs() {return cl_probs;}
	int                     get_layer_amount(){ return layer_amount; }
	shared_ptr <int[]>      get_init_layer_idx (){return init_layer_idx;}
	shared_ptr <int[]>      get_init_layer_size() { return layer_size; }
	~initial_prob_img() {
		
		delete []init_prob_img;
	};
};

