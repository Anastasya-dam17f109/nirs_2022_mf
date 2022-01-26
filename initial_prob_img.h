#pragma once
#include "mix_img_obj.h"


class initial_prob_img
{
	shared_ptr<mix_img_obj> m_image;
	shared_ptr <int[]>      layer_size;
	int layer_amount = 1;
	shared_ptr <shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>[]> init_prob_img;
	shared_ptr<double[]> mix_shift;
	shared_ptr<double[]> mix_scale;
	shared_ptr<double[]> cl_probs;
	unsigned image_len_x = 32;
	unsigned image_len_y = 32;
	int class_amount = 1;

public:
	initial_prob_img(int img_size, mix_type mix_t, int amount_targets, int classes);
	initial_prob_img(shared_ptr<mix_img_obj> image);
	initial_prob_img(string filename);
	void generate_init_probs_mixtures();
	int      get_class_amount() { return class_amount; }
	std::pair<int, int>   get_image_len() { return std::pair<int, int>(image_len_x, image_len_y); }
	shared_ptr <shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>[]>  get_image() { return init_prob_img;}
	shared_ptr<double[]> get_cl_probs() {return cl_probs;}
	~initial_prob_img() {};
};

