//#include "pch.h"
#include "initial_prob_img.h"


initial_prob_img::initial_prob_img(shared_ptr<mix_img_obj> image){
	m_image      = image;
	image_len_x  = m_image->get_image_len().first;
	image_len_y  = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift    = m_image->get_shift();
	mix_scale    = m_image->get_scale();
	cl_probs     = shared_ptr<double[]>(new double[class_amount]);
	layer_amount = m_image->get_layer_amount();
	layer_size   = m_image->get_layer_size();
	layer_idx    = m_image->get_layer_idx();
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;
	alloc_layer_mmr();
	generate_init_probs_mixtures();
}

// создание объекта через загрузчик-файл

initial_prob_img::initial_prob_img(string filename) {
	m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(filename));
	image_len_x  = m_image->get_image_len().first;
	image_len_y  = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift    = m_image->get_shift();
	mix_scale    = m_image->get_scale();
	layer_amount = m_image->get_layer_amount();
	layer_size   = m_image->get_layer_size();
	layer_idx = m_image->get_layer_idx();
	cl_probs     = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;
	alloc_layer_mmr();
	generate_init_probs_mixtures();
}

// создание объекта через заданные значения для генератора картинки

initial_prob_img::initial_prob_img(int img_size, mix_type mix_t, int amount_targets, int classes) {
	m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(img_size, mix_t,amount_targets, classes));
	image_len_x  = m_image->get_image_len().first;
	image_len_y  = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift    = m_image->get_shift();
	mix_scale    = m_image->get_scale();
	layer_amount = m_image->get_layer_amount();
	layer_size   = m_image->get_layer_size();
	layer_idx = m_image->get_layer_idx();;
	cl_probs     = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;

	alloc_layer_mmr();
	generate_init_probs_mixtures();
}

//

initial_prob_img::initial_prob_img(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes) {
	m_image = shared_ptr<mix_img_obj>(new mix_img_obj(file_name, img_size, mix_t, amount_targets, classes));
	image_len_x = m_image->get_image_len().first;
	image_len_y = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift = m_image->get_shift();
	mix_scale = m_image->get_scale();
	layer_amount = m_image->get_layer_amount();
	layer_size = m_image->get_layer_size();
	cl_probs = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;

	alloc_layer_mmr();
	generate_init_probs_mixtures();
}

void initial_prob_img::alloc_layer_mmr(){

	init_layer_idx = shared_ptr <int[]>(new int [layer_amount]);
	int summ = 0;
	for (int k = 0; k < layer_amount; ++k) {
		init_layer_idx[k] = summ;
		summ += layer_size[k] * layer_size[k] * class_amount;
	}

	init_prob_img = new double [summ];
	
}

// вычисление начальных значений вероятностей через одномерные смеси

void initial_prob_img::generate_init_probs_mixtures() {
	double pi = 3.14;
	double summ1;
	for (int k = 0; k < layer_amount; ++k) {
		for (int i = 0; i < layer_size[k]; i++) {
			for (int j = 0; j < layer_size[k]; j++) {
				summ1 = 0;
				for (int t = 0; t < class_amount; t++) {
					//#pragma omp critical
					//{
					if (mix_scale[t] != 0)
						if (m_image->get_mixture_type() == NORMAL)
							summ1 += (1 / (mix_scale[t] * sqrt(2 * pi)))*exp(-(pow(m_image->get_image()[layer_idx[k] +
								i * layer_size[k] + 
								j] 
								- mix_shift[t], 2)) /
									(2.0 * mix_scale[t] * mix_scale[t]));
						else {
							if (m_image->get_mixture_type() == LOGNORMAL)
								summ1 += (1 / (mix_scale[t] * m_image->get_image()[layer_idx[k]+i* layer_size[k] +j] * sqrt(2 * pi)))
												*exp(-(pow(log(m_image->get_image()[layer_idx[k] + i * layer_size[k] + j]) - mix_shift[t], 2)) /
														(2.0 * mix_scale[t] * mix_scale[t]));
							else {
								if (m_image->get_mixture_type() == RAYLEIGH)
								summ1 += (1 / (mix_scale[t] * mix_scale[t]))*exp(-(pow(m_image->get_image()[layer_idx[k] + i * layer_size[k] + j], 2)) /
									(2.0 * mix_scale[t] * mix_scale[t]));
							}
						}
					//}
				}

				for (int t = 0; t < class_amount; t++) {
					if (mix_scale[t] != 0)
						if (m_image->get_mixture_type() == NORMAL)
							init_prob_img[init_layer_idx[k] + i* layer_size[k]*  class_amount +j * class_amount +t] =
								(1 / (mix_scale[t] * sqrt(2 * pi)*summ1))
								*exp(-(pow(m_image->get_image()[layer_idx[k] + i * layer_size[k] + j] - mix_shift[t], 2))
								/ (2.0 * mix_scale[t] * mix_scale[t]));
						else {
							if (m_image->get_mixture_type() == LOGNORMAL)
								init_prob_img[init_layer_idx[k] + 
								i * layer_size[k] * class_amount + 
								j * class_amount + 
								t]
									= (1 / (mix_scale[j] * m_image->get_image()[layer_idx[k] + i * layer_size[k] + j] * sqrt(2 * pi)*summ1))*
											exp(-(pow(log(m_image->get_image()[layer_idx[k] + i * layer_size[k] + j]) - mix_shift[t], 2))
													/ (2.0 * mix_scale[j] * mix_scale[t]));
							else {
								if (m_image->get_mixture_type() == RAYLEIGH)
									init_prob_img[init_layer_idx[k] +
									i * layer_size[k] * class_amount + 
									j * class_amount +
									t]
									= (1 / (mix_scale[t] * mix_scale[t] *summ1))
										*exp(-(pow(m_image->get_image()[layer_idx[k] + i * layer_size[k] + j], 2))
										/ (2.0 * mix_scale[t] * mix_scale[t]));
							}
						}
				}
			}
		}
	}
}



