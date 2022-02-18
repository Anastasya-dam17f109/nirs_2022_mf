#pragma once
#include "mix_img_obj.h"
#include "omp.h"
#include "quad_tree_handler.h"

//класс-обработчик большого изображения

class mixture_handler
{
	shared_ptr<initial_prob_img> m_image;
	double** raw_image;
	unsigned** class_flag;
	shared_ptr<string[]> img_mask_list;

	unsigned img_l_x;
	unsigned img_l_y;

	shared_ptr<double[]> mix_shift;
	shared_ptr<double[]> mix_scale;
	shared_ptr<double[]> mix_weight;
	shared_ptr<double[]> mix_prob;

	unsigned window_size;
	unsigned min_trg_size = 17;

	double accuracy = 0.001;

	unsigned hyp_cl_amount = 1;
	unsigned hyp_cl_amount_mod = 1;
	bool equal_hyp_flag = true;

	shared_ptr<shared_ptr<double[]>[]> g_i_j;
	double bic_value = 0;
	shared_ptr<double[]> rfar;
	float all_mistakes = 0;
	shared_ptr<double[]> mistake_mix;
	mix_type mixture_type;
	std::ofstream out;
	string gen_mix_filename = "D:\\generated_image.txt";
	string split_mix_filename = "D:\\splitted_image.txt";

public:
	mixture_handler(shared_ptr < mix_img_obj> img, unsigned h_classes, double acc);
	mixture_handler(shared_ptr <initial_prob_img> img, unsigned h_classes, double acc);
	//mixture_handler(shared_ptr<shared_ptr<double[]>[]> img, unsigned h_classes, double acc);
	void draw_graphics();
	void detect_result_by_mask();
	
	bool get_equal_hyp_flag() {
		return equal_hyp_flag;
	}
	double get_bic_value() {return bic_value;}
	void printInformation();
	void printInformation_to_image();
	void mixture_inicalization();
	
	double find_k_stat(double * data, int wind_size, int k_stat);
	double find_med(double* window, int wind_size);
	std::pair<int, int> partition(double* mass, int left, int right, int  ind_pivot);
	void  quickSort(double * data, int wind_size, int l, int r);

	
    void mixture_optimal_redraw_opMP();
    void mixture_optimal_redraw_opMP_V2();
	void kolmogorov_optimal_redraw_opMP();
    void q_tree_optimal_redraw_opMP();
	
	void BIC();

	~mixture_handler() {
        for (unsigned i = 0; i < img_l_x; i++)
			delete[] class_flag[i];
		delete[] class_flag;
	}
};

