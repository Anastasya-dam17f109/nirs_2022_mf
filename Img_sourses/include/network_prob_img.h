#pragma once
#include "mix_img_obj.h"
#include "omp.h"

class network_prob_img
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

	int init_window_size = 8;

	bool genFlag = true;
public:
	// конструктор получает имя файла конфигурации - имена всех файлов, что сделала нейросеть
	network_prob_img(string configFilename);
	void     alloc_layer_mmr();

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

