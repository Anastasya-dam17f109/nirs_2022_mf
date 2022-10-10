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
	//generate_init_probs_mixtures();
	generate_init_probs_mix_opMP_V2();

}

// создание объекта через загрузчик-файл +  флаг, необходимо ли считать статистики по изображениям и генерить вероятности

initial_prob_img::initial_prob_img(string filename, bool flag) {
	genFlag = flag;
	m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(filename, genFlag));

	image_len_x  = m_image->get_image_len().first;
	image_len_y  = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();

	if (genFlag)
	{
		mix_shift = m_image->get_shift();
		mix_scale = m_image->get_scale();
	}
	layer_amount = m_image->get_layer_amount();
	layer_size   = m_image->get_layer_size();
	layer_idx    = m_image->get_layer_idx();
	cl_probs     = shared_ptr<double[]>(new double[class_amount]);
	for (int i = 0; i < class_amount; ++i)
		cl_probs[i] = 1.0 / class_amount;
	alloc_layer_mmr();
	if (genFlag)
	{
		//generate_init_probs_mixtures();
		generate_init_probs_mix_opMP_V2();
		//mix_scale = m_image->get_scale();
	}
}

//

initial_prob_img::initial_prob_img(string configFilename)
{
	

	/*m_image = shared_ptr<mix_img_obj>(new mix_img_obj(filename_load_image, false));

	image_len_x = m_image->get_image_len().first;
	image_len_y = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();*/
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
	//generate_init_probs_mix_opMP_V2();
	//generate_init_probs_max_apost_opMP_V2();
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
	//generate_init_probs_mixtures();
	//generate_init_probs_mix_opMP_V2();
	generate_init_probs_max_apost_opMP_V2();
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

// вычисление начальных значений вероятностей через обработку EM алгоритмом  по локальным областям

void initial_prob_img::generate_init_probs_mix_opMP_V2() {

	int window_size = init_window_size;
	float accuracy = 0.001;
	//for (int i = 0; i < class_amount; ++i)
		//mix_scale[i] = mix_scale[i] / (pow(2, layer_amount));
	for (int k = layer_amount-1; k >-1; --k) {
		window_size = init_window_size;
		int add_amount_x = layer_size[k] % window_size;
		int add_amount_y = layer_size[k] % window_size;
		int amount_window_x = layer_size[k] / window_size;
		int amount_window_y = layer_size[k] / window_size;
		if (amount_window_x == 0){
			amount_window_x = 1;
			amount_window_y = 1;
			add_amount_x = 0;
			add_amount_y = 0;
			window_size = layer_size[k];
		}
		int u_new_n = (window_size + add_amount_x) *(window_size + add_amount_y);
		int thr_nmb = 4;
		double** new_g_ij = new double*[u_new_n *thr_nmb];
		double** new_g_ij_0 = new double*[u_new_n * thr_nmb];
		
		for (int l = 0; l < u_new_n * thr_nmb; ++l) {
			new_g_ij[l] = new double[class_amount];
			new_g_ij_0[l] = new double[class_amount];
			for (int t = 0; t < class_amount; t++) {
				new_g_ij[l][t] = 0;
				new_g_ij_0[l][t] = 0;
			}
		}
		
		
		auto begin1 = std::chrono::steady_clock::now();
#pragma omp parallel
		{
			int loc_window_size = window_size;
			int loc_hyp_cl_amount = class_amount;
			int x_l = loc_window_size;
			int y_l = loc_window_size;
			int ofset = omp_get_thread_num();
			int itr, x_min, y_min, j, loc_u_new_n, t, l;
			//x_min = init_layer_idx[k]; y_min = init_layer_idx[k];
			double * new_weights = new  double[loc_hyp_cl_amount];
			double * buf_new_weights = new  double[loc_hyp_cl_amount];

			for (l = 0; l < loc_hyp_cl_amount; ++l) 
				new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
			
			
			double summ = 0;
			const double sq_pi = sqrt(2 * pi);
			double pix_buf, cur_max;
			double last_cur_max = 0;

			double buf_max = 0;
			bool stop_flag = true;
			unsigned idx_max = 0;
#pragma omp for
			for (int r = 0; r < amount_window_x; ++r) {
				x_min = r * loc_window_size;
				if (r < amount_window_x - 1)
					x_l = loc_window_size;
				else
					x_l = loc_window_size + add_amount_x;
				for (j = 0; j < amount_window_y; ++j) {

					y_min = j * loc_window_size;
					if (j < amount_window_y - 1)
						y_l = loc_window_size;
					else
						y_l = loc_window_size + add_amount_y;


					itr = 0;
					stop_flag = true;
					cur_max = 0;
					loc_u_new_n = y_l * x_l;

					while (stop_flag && (itr < 500)) {
						++itr;

						for (l = ofset * u_new_n; l < ofset* u_new_n + loc_u_new_n; ++l) {
							summ = 0;
							
							pix_buf = m_image->get_image()[layer_idx[k] +
								(x_min + (l - ofset * u_new_n) / y_l)*layer_size[k] + y_min + (l - ofset * u_new_n) % y_l];
							for (t = 0; t < loc_hyp_cl_amount; ++t) {
								if (m_image->get_mixture_type() == NORMAL)
									summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi))*exp(-((pix_buf
										- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
								else {
									if (m_image->get_mixture_type() == RAYLEIGH)
										summ += new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-((pix_buf
											)*(pix_buf)) / (2.0 * mix_scale[t] * mix_scale[t]));

									else {
										if (m_image->get_mixture_type() == LOGNORMAL)
											summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi*pix_buf))*exp(-((log(pix_buf)
												- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
									}
								}
								/*cout << "l  " << summ << " " << pix_buf << " " << layer_idx[k] +
									(x_min + (l - ofset * u_new_n) / y_l)*layer_size[k] + y_min + (l - ofset * u_new_n) % y_l 
									<< " " << layer_idx[layer_amount-1]+1024*1024 << endl;*/

							}

							for (t = 0; t < loc_hyp_cl_amount; ++t) {
								if (l == ofset * u_new_n)
									buf_new_weights[t] = 0;
								if (m_image->get_mixture_type() == NORMAL) {
									
									new_g_ij[l][t] = new_weights[t] * (1 / (mix_scale[t] * sq_pi*summ))*exp(-((pix_buf
										- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
									/*cout << "l  " << new_g_ij[l][t] << endl;*/

									//if (layer_size[k] == 32)
									//	cout << "l  " << l << " " << t << " " << (1 / (mix_scale[t] * sq_pi*summ))*exp(-((pix_buf
									//		- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t])) << " " << exp(-((pix_buf
									//		- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t])) << " " << summ <<   endl;
								}

								else {
									if (m_image->get_mixture_type() == RAYLEIGH)
										new_g_ij[l][t] = new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-(pix_buf
											*pix_buf) / (2.0 * mix_scale[t] * mix_scale[t]));

									else {
										if (m_image->get_mixture_type() == LOGNORMAL) {
											//cout << "l  " << l << endl;
											new_g_ij[l][t] = new_weights[t] * (1 / (mix_scale[t] * sq_pi*summ*pix_buf))*exp(-((log(pix_buf)
												- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));

										}
									}
								}
								buf_new_weights[t] += new_g_ij[l][t];
								if (l == ofset * u_new_n + loc_u_new_n - 1)
									new_weights[t] = buf_new_weights[t] / double(loc_u_new_n);
								if (cur_max < abs(new_g_ij[l][t] - new_g_ij_0[l][t]))
									cur_max = abs(new_g_ij[l][t] - new_g_ij_0[l][t]);
								new_g_ij_0[l][t] = new_g_ij[l][t];
							}

						}

						if (stop_flag) {
							if (cur_max != 0)
								last_cur_max = cur_max;
							(cur_max < accuracy) ? stop_flag = false : cur_max = 0;
						}
					}
#pragma omp critical
					{
						//if (layer_size[k] == 32)
						//cout << "jjjjj" << endl;

						for (t = ofset * u_new_n; t < ofset * u_new_n + loc_u_new_n; ++t) {
							
							for (l = 0; l < loc_hyp_cl_amount; ++l) {
								if (t == ofset * u_new_n)
									//new_weights[l] = mix_prob[l];
									new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
								init_prob_img[init_layer_idx[k] +
									(x_min + (t - ofset * u_new_n) / y_l) * layer_size[k] * class_amount +
									(y_min + (t - ofset * u_new_n) % y_l) * class_amount +
									l] = new_g_ij[t][l];
								//if (layer_size[k] == 32)
									//cout << t << " " << l << " " << new_g_ij[t][l] << endl;
								new_g_ij_0[t][l] = 0;
							}
							
						}
					}
				}
			}
			delete[] new_weights;
			delete[] buf_new_weights;
			
		}
		/*for (int i = 0; i < class_amount; ++i)
			mix_shift[i] = mix_shift[i] * 2;*/
		for (int t = 0; t < u_new_n* thr_nmb; ++t) {
			delete[] new_g_ij[t];
			delete[] new_g_ij_0[t];
		}

		delete[] new_g_ij;
		delete[] new_g_ij_0;
		auto end1 = std::chrono::steady_clock::now();
		auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
		cout << "elapsed_ms1  " << elapsed_ms1.count() << endl;

	}
}


void initial_prob_img::generate_init_probs_max_apost_opMP_V2() {
	int window_size = init_window_size;
	float accuracy = 0.001;
	for (int i = 0; i < class_amount; ++i)
		mix_scale[i] = mix_scale[i] / (pow(2, layer_amount));

	for (int k = 0; k < layer_amount; ++k) {
		window_size = init_window_size;
		int add_amount_x = layer_size[k] % window_size;
		int add_amount_y = layer_size[k] % window_size;
		int amount_window_x = layer_size[k] / window_size;
		int amount_window_y = layer_size[k] / window_size;
		if (amount_window_x == 0) {
			amount_window_x = 1;
			amount_window_y = 1;
			add_amount_x = 0;
			add_amount_y = 0;
			window_size = layer_size[k];
		}
		int u_new_n = (window_size + add_amount_x) *(window_size + add_amount_y);
		int thr_nmb = 4;
		double* buf_mass = new double[u_new_n *thr_nmb];
		
		for (int i = 0; i < class_amount; ++i)
			mix_scale[i] = mix_scale[i] *2;

		auto begin1 = std::chrono::steady_clock::now();
#pragma omp parallel
		{
			auto L_max_calculation = [&](double* data, int data_size, int beg, 
										int iter, double _mix_shift, double _mix_scale, double* max_L_mass) {
				double buf_max_l = 0;
				bool flag = false;
				double B;
				for (int m = beg; m < beg + data_size; m++) {
					B = 0;
					if (mix_scale != 0) {
						if (m_image->get_mixture_type() == NORMAL)
							B = (1.0 / (_mix_scale))*
							exp(-(pow(data[m] - _mix_shift, 2)) /
							(2.0 * _mix_scale * _mix_scale));
						else
							if (m_image->get_mixture_type() == LOGNORMAL)
								B = (1.0 / (_mix_scale*data[m]))*
								exp(-(pow(log(data[m]) - _mix_shift, 2)) /
								(2.0 * _mix_scale * _mix_scale));
							else
								if (m_image->get_mixture_type() == RAYLEIGH)
									B = (data[m] / pow(_mix_scale, 2))*
									exp(-(pow(data[m], 2)) /
									(2.0 * _mix_scale * _mix_scale));
					}
					else {
						flag = true;
						break;
					}

					if (B > 0)
						buf_max_l = buf_max_l + log(B / (data_size));
				}

				max_L_mass[iter] = buf_max_l;
			};
			int loc_window_size = window_size;
			int loc_hyp_cl_amount = class_amount;
			int x_l = loc_window_size;
			int y_l = loc_window_size;
			int ofset = omp_get_thread_num();
			int itr, x_min, y_min, j, loc_u_new_n, t, l;
			int max_idx;
			double * max_L_mass = new  double[loc_hyp_cl_amount];

			double pix_buf, cur_max;
			

			double buf_max = 0;
			bool stop_flag = true;
			unsigned idx_max = 0;
			#pragma omp for
			for (int r = 0; r < amount_window_x; ++r) {
				x_min = r * loc_window_size;
				if (r < amount_window_x - 1)
					x_l = loc_window_size;
				else
					x_l = loc_window_size + add_amount_x;
				for (j = 0; j < amount_window_y; ++j) {

					y_min = j * loc_window_size;
					if (j < amount_window_y - 1)
						y_l = loc_window_size;
					else
						y_l = loc_window_size + add_amount_y;
					loc_u_new_n = y_l * x_l;

					for (l = ofset * u_new_n; l < ofset* u_new_n + loc_u_new_n; ++l) {
						buf_mass[l] = m_image->get_image()[layer_idx[k] +
							(x_min + (l - ofset * u_new_n) / y_l)*layer_size[k] + y_min + (l - ofset * u_new_n) % y_l];
					}
					for (t = 0; t < loc_hyp_cl_amount; ++t) {
						L_max_calculation(buf_mass, loc_u_new_n, ofset * u_new_n,
							t, mix_shift[t], mix_scale[t], max_L_mass);
						if (t == 0) {
							cur_max = max_L_mass[0];
							max_idx = 0;
						}
						else {
							if (cur_max < max_L_mass[t]) {
								cur_max = max_L_mass[t];
								max_idx = t;
							}
						}

					}
					#pragma omp critical
					{
						for (t = ofset * u_new_n; t < ofset * u_new_n + loc_u_new_n; ++t) {

							for (l = 0; l < loc_hyp_cl_amount; ++l) {
								
								init_prob_img[init_layer_idx[k] +
									(x_min + (t - ofset * u_new_n) / y_l) * layer_size[k] * class_amount +
									(y_min + (t - ofset * u_new_n) % y_l) * class_amount +
									l] = max_idx + 1;
								
							}

						}
					}
				}
			}
			delete[] max_L_mass;
			

		}

		delete[] buf_mass;
		auto end1 = std::chrono::steady_clock::now();
		auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
		cout << "elapsed_ms1  " << elapsed_ms1.count() << endl;
	}
}