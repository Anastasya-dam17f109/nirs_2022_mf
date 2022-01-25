#include "pch.h"
#include "initial_prob_img.h"


initial_prob_img::initial_prob_img(shared_ptr<mix_img_obj> image){
	m_image = image;
	image_len_x  = m_image->get_image_len().first;
	image_len_y  = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift = m_image->get_shift();
	mix_scale = m_image->get_scale();
	cl_probs = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;

	init_prob_img = shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>(new shared_ptr<shared_ptr<double[]>[]>[image_len_x]);
	for (int i = 0; i < image_len_x; i++) {
		init_prob_img[i] = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[image_len_y]);
		for (int j = 0; j < image_len_y; j++)
			init_prob_img[i][j] = shared_ptr<double[]>(new double[class_amount]);
	}
	generate_init_probs_mixtures();
}

// создание объекта через загрузчик-файл

initial_prob_img::initial_prob_img(string filename) {
	//mix_img_obj *img = new mix_img_obj(filename);

	m_image = shared_ptr<mix_img_obj>(new mix_img_obj(filename));
	image_len_x = m_image->get_image_len().first;
	image_len_y = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift = m_image->get_shift();
	mix_scale = m_image->get_scale();
	cl_probs = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;

	init_prob_img = shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>(new shared_ptr<shared_ptr<double[]>[]>[image_len_x]);
	for (int i = 0; i < image_len_x; i++) {
		init_prob_img[i] = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[image_len_y]);
		for (int j = 0; j < image_len_y; j++)
			init_prob_img[i][j] = shared_ptr<double[]>(new double[class_amount]);
	}
	generate_init_probs_mixtures();
}

// создание объекта через заданные значения для генератора картинки

initial_prob_img::initial_prob_img(int img_size, mix_type mix_t, int amount_targets, int classes) {
	m_image = shared_ptr<mix_img_obj>(new mix_img_obj(img_size, mix_t,amount_targets, classes));
	image_len_x = m_image->get_image_len().first;
	image_len_y = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
	mix_shift = m_image->get_shift();
	mix_scale = m_image->get_scale();
	cl_probs = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;


	init_prob_img = shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>(new shared_ptr<shared_ptr<double[]>[]>[image_len_x]);
	for (int i = 0; i < image_len_x; i++) {
		init_prob_img[i] = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[image_len_y]);
		for (int j = 0; j < image_len_y; j++)
			init_prob_img[i][j] = shared_ptr<double[]>(new double[class_amount]);
	}
	generate_init_probs_mixtures();
}

// вычисление начальных значений вероятностей через одномерные смеси

void initial_prob_img::generate_init_probs_mixtures() {
	double pi = 3.14;
	double summ1;
	for (int i = 0; i < image_len_x; i++) {
		for (int j = 0; j < image_len_y; j++) {
			summ1 = 0;
			for (int t = 0; t < class_amount; t++) {
				//#pragma omp critical
				//{
					if (mix_scale[t] != 0)
						if (m_image->get_mixture_type() == NORMAL)
							summ1 += (1 / (mix_scale[t] * sqrt(2 * pi)))*exp(-(pow(m_image->get_image()[i][j] - mix_shift[t], 2)) /
							(2.0 * mix_scale[t] * mix_scale[t]));
						else {
							if (m_image->get_mixture_type() == LOGNORMAL)
								summ1 += (1 / (mix_scale[t] * m_image->get_image()[i][j] * sqrt(2 * pi)))*exp(-(pow(log(m_image->get_image()[i][j]) - mix_shift[t], 2)) /
								(2.0 * mix_scale[t] * mix_scale[t]));
						}
				//}
			}

			for (int t = 0; t < class_amount; t++) {
				if (mix_scale[t] != 0)
					if (m_image->get_mixture_type() == NORMAL)
						init_prob_img[i][j][t] = (1 / (mix_scale[t] * sqrt(2 * pi)*summ1))*exp(-(pow(m_image->get_image()[i][j] - mix_shift[t], 2))
							/ (2.0 * mix_scale[t] * mix_scale[t]));
					else {
						if (m_image->get_mixture_type() == LOGNORMAL)
							init_prob_img[i][j][t] = (1 / (mix_scale[j] * m_image->get_image()[i][j] * sqrt(2 * pi)*summ1))*exp(-(pow(log(m_image->get_image()[i][j]) - mix_shift[t], 2))
								/ (2.0 * mix_scale[j] * mix_scale[t]));
					}
			}
		}
	}
}



