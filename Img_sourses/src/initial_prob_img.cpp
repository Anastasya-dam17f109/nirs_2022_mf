#include "initial_prob_img.h"


void initial_prob_img::gen_prob_img_from_obj(shared_ptr<mix_img_obj> image){
    m_image      = image;
    m_image_len_x  = m_image->get_image_len().first;
    m_image_len_y  = m_image->get_image_len().second;
    m_class_amount = m_image->get_class_amount();
    mix_shift    = m_image->get_shift();
    mix_scale    = m_image->get_scale();
    m_cl_probs     = shared_ptr<double[]>(new double[m_class_amount]);
    m_layer_amount = m_image->get_layer_amount();
    m_layer_size   = m_image->get_layer_size();
    m_layer_idx    = m_image->get_layer_idx();
    for (int i = 0; i < m_class_amount; ++i)
        m_cl_probs[i] = 1.0 / m_class_amount;
    alloc_layer_mmr();
    //generate_init_probs_mixtures();
    generate_init_probs_mix_opMP_V2();

}

// создание объекта через загрузчик-файл +  флаг, необходимо ли считать статистики по изображениям и генерить вероятности

void initial_prob_img::gen_prob_img_from_config(string filename, int mode) {
    m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(filename, true, mode));
    m_image_len_x  = m_image->get_image_len().first;
    m_image_len_y  = m_image->get_image_len().second;
    m_class_amount = m_image->get_class_amount();
    mix_shift = m_image->get_shift();
    mix_scale = m_image->get_scale();
    m_layer_amount = m_image->get_layer_amount();
    m_layer_size   = m_image->get_layer_size();
    m_layer_idx    = m_image->get_layer_idx();
    m_cl_probs     = shared_ptr<double[]>(new double[m_class_amount]);
    for (int i = 0; i < m_class_amount; ++i)
        m_cl_probs[i] = 1.0 / m_class_amount;
    alloc_layer_mmr();
    generate_init_probs_mix_opMP_V2();
}


// создание объекта через заданные значения для генератора картинки

void initial_prob_img::gen_prob_img_with_targs(int img_size, mix_type mix_t, int amount_targets, int classes) {
    m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(img_size, mix_t,amount_targets, classes));
    m_image_len_x  = m_image->get_image_len().first;
    m_image_len_y  = m_image->get_image_len().second;
    m_class_amount = m_image->get_class_amount();
    mix_shift    = m_image->get_shift();
    mix_scale    = m_image->get_scale();
    m_layer_amount = m_image->get_layer_amount();
    m_layer_size   = m_image->get_layer_size();
    m_layer_idx = m_image->get_layer_idx();;
    m_cl_probs     = shared_ptr<double[]>(new double[m_class_amount]);
    for (int i = 0; i < m_class_amount; ++i)
        m_cl_probs[i] = 1.0 / m_class_amount;

    alloc_layer_mmr();
    generate_init_probs_mixtures();
}

//

void initial_prob_img::gen_prob_img_with_targs_from_file(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes) {
    m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(file_name, img_size, mix_t, amount_targets, classes));
    m_image_len_x  = m_image->get_image_len().first;
    m_image_len_y  = m_image->get_image_len().second;
    m_class_amount = m_image->get_class_amount();
    mix_shift    = m_image->get_shift();
    mix_scale    = m_image->get_scale();
    m_layer_amount = m_image->get_layer_amount();
    m_layer_size = m_image->get_layer_size();
    m_cl_probs = shared_ptr<double[]>(new double[m_class_amount]);
    for (int i = 0; i < m_class_amount; ++i)
        m_cl_probs[i] = 1.0 / m_class_amount;

    alloc_layer_mmr();

    generate_init_probs_max_apost_opMP_V2();
}


// вычисление начальных значений вероятностей через одномерные смеси

void initial_prob_img::generate_init_probs_mixtures() {
    double pi = 3.14;
    double summ1;
    for (int k = 0; k < m_layer_amount; ++k) {
        for (int i = 0; i < m_layer_size[2*k]; i++) {
            for (int j = 0; j < m_layer_size[2*k +1]; j++) {
                summ1 = 0;
                for (int t = 0; t < m_class_amount; t++) {
                    //#pragma omp critical
                    //{
                    if (mix_scale[t] != 0)
                        if (m_image->get_mixture_type() == NORMAL)
                            summ1 += (1 / (mix_scale[t] * sqrt(2 * pi)))*exp(-(pow(m_image->get_image()[m_layer_idx[k] +
                                                                               i * m_layer_size[2 * k + 1] +
                                                                               j]
                                                                               - mix_shift[t], 2)) /
                                    (2.0 * mix_scale[t] * mix_scale[t]));
                        else {
                            if (m_image->get_mixture_type() == LOGNORMAL)
                                summ1 += (1 / (mix_scale[t] * m_image->get_image()[m_layer_idx[k]+i* m_layer_size[2 * k + 1] +j] * sqrt(2 * pi)))
                                        *exp(-(pow(log(m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j]) - mix_shift[t], 2)) /
                                        (2.0 * mix_scale[t] * mix_scale[t]));
                            else {
                                if (m_image->get_mixture_type() == RAYLEIGH)
                                    summ1 += (1 / (mix_scale[t] * mix_scale[t]))*exp(-(pow(m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j], 2)) /
                                            (2.0 * mix_scale[t] * mix_scale[t]));
                            }
                        }
                    //}
                }

                for (int t = 0; t < m_class_amount; t++) {
                    if (mix_scale[t] != 0)
                        if (m_image->get_mixture_type() == NORMAL)
                            m_prob_img[m_init_layer_idx[k] + i* m_layer_size[2 * k + 1]*  m_class_amount +j * m_class_amount +t] =
                                    (1 / (mix_scale[t] * sqrt(2 * pi)*summ1))
                                    *exp(-(pow(m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j] - mix_shift[t], 2))
                                    / (2.0 * mix_scale[t] * mix_scale[t]));
                        else {
                            if (m_image->get_mixture_type() == LOGNORMAL)
                                m_prob_img[m_init_layer_idx[k] +
                                        i * m_layer_size[2 * k + 1] * m_class_amount +
                                        j * m_class_amount +
                                        t]
                                        = (1 / (mix_scale[j] * m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j] * sqrt(2 * pi)*summ1))*
                                        exp(-(pow(log(m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j]) - mix_shift[t], 2))
                                        / (2.0 * mix_scale[j] * mix_scale[t]));
                            else {
                                if (m_image->get_mixture_type() == RAYLEIGH)
                                    m_prob_img[m_init_layer_idx[k] +
                                            i * m_layer_size[2 * k + 1] * m_class_amount +
                                            j * m_class_amount +
                                            t]
                                            = (1 / (mix_scale[t] * mix_scale[t] *summ1))
                                            *exp(-(pow(m_image->get_image()[m_layer_idx[k] + i * m_layer_size[2 * k + 1] + j], 2))
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

    int window_size = m_init_window_size;
    float accuracy = 0.001;
	const double sq_pi = sqrt(2 * pi);
    
    for (int k = m_layer_amount-1; k >-1; --k) {
        window_size = m_init_window_size;
        int add_amount_x = m_layer_size[2*k] % window_size;
        int add_amount_y = m_layer_size[2 * k + 1] % window_size;
        int amount_window_x = m_layer_size[2*k] / window_size;
        int amount_window_y = m_layer_size[2 * k + 1] / window_size;
        if (amount_window_x == 0){
            amount_window_x = 1;
            amount_window_y = 1;
            add_amount_x = 0;
            add_amount_y = 0;
            window_size = m_layer_size[2*k];
        }
        int u_new_n = (window_size + add_amount_x) * (window_size + add_amount_y);
        int thr_nmb = 4;
        double** new_g_ij = new double * [u_new_n * thr_nmb];
        double** new_g_ij_0 = new double * [u_new_n * thr_nmb];
		double * new_weights = new double[m_class_amount * thr_nmb];
		double * buf_new_weights = new double[m_class_amount * thr_nmb];

        for (int l = 0; l < u_new_n * thr_nmb; ++l) {
            new_g_ij[l] = new double[m_class_amount];
            new_g_ij_0[l] = new double[m_class_amount];
            for (int t = 0; t < m_class_amount; t++) {
                new_g_ij[l][t] = 0;
                new_g_ij_0[l][t] = 0;
            }
        }

        auto begin1 = std::chrono::steady_clock::now();
		//if ( (m_layer_size[2 * k] > 4))
		{
#pragma omp parallel
			{
				int x_l = window_size, y_l = window_size, ofset = omp_get_thread_num();
				int itr, x_min, y_min, j, loc_u_new_n, t, l;
				double pix_buf, cur_max, summ = 0, last_cur_max = 0, buf_max = 0;
				bool stop_flag = true;
				unsigned idx_max = 0;
				for (l = ofset * m_class_amount; l < ofset * m_class_amount + m_class_amount; ++l)
					new_weights[l] = 1.0 / double(m_class_amount);
#pragma omp for
				for (int r = 0; r < amount_window_x; ++r) {
					x_min = r * window_size;
					if (r < amount_window_x - 1)
						x_l = window_size;
					else
						x_l = window_size + add_amount_x;

					for (j = 0; j < amount_window_y; ++j) {
						y_min = j * window_size;

						if (j < amount_window_y - 1)
							y_l = window_size;
						else
							y_l = window_size + add_amount_y;

						itr = 0;
						stop_flag = true;
						cur_max = 0;
						loc_u_new_n = y_l * x_l;

						while (stop_flag && (itr < 500)) {
							++itr;

							for (l = ofset * u_new_n; l < ofset * u_new_n + loc_u_new_n; ++l) {
								summ = 0;

								pix_buf = m_image->get_image()[m_layer_idx[k] +
									(x_min + (l - ofset * u_new_n) / y_l)*m_layer_size[2 * k + 1] + y_min + (l - ofset * u_new_n) % y_l];
								for (t = 0; t < m_class_amount; ++t) {
									if (m_image->get_mixture_type() == NORMAL)
										summ += new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi))
										*exp(-((pix_buf - mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
									else {
										if (m_image->get_mixture_type() == RAYLEIGH)
											summ += new_weights[ofset * m_class_amount + t] * (pix_buf / (mix_scale[t] * mix_scale[t]))
											* exp(-((pix_buf)*(pix_buf)) / (2.0 * mix_scale[t] * mix_scale[t]));

										else {
											if (m_image->get_mixture_type() == LOGNORMAL)
												summ += new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*pix_buf))*exp(-((log(pix_buf)
													- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
										}
									}
									/*cout << "l  " << summ << " " << pix_buf << " " << layer_idx[k] +
										(x_min + (l - ofset * u_new_n) / y_l)*layer_size[k] + y_min + (l - ofset * u_new_n) % y_l
										<< " " << layer_idx[layer_amount-1]+1024*1024 << endl;*/

								}

								for (t = 0; t < m_class_amount; ++t) {
									if (l == ofset * u_new_n)
										buf_new_weights[ofset * m_class_amount + t] = 0;
									if (m_image->get_mixture_type() == NORMAL)
										new_g_ij[l][t] = new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*summ))
										* exp(-((pix_buf - mix_shift[t]) * (pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));

									else {
										if (m_image->get_mixture_type() == RAYLEIGH)
											new_g_ij[l][t] = new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))
											* exp(-(pix_buf * pix_buf) / (2.0 * mix_scale[t] * mix_scale[t]));

										else {
											if (m_image->get_mixture_type() == LOGNORMAL) {
												//cout << "l  " << l << endl;
												new_g_ij[l][t] = new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*summ*pix_buf))
													* exp(-((log(pix_buf) - mix_shift[t]) * (log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));

											}
										}
									}
									buf_new_weights[ofset * m_class_amount + t] += new_g_ij[l][t];
									if (l == ofset * u_new_n + loc_u_new_n - 1)
										new_weights[ofset * m_class_amount + t] = buf_new_weights[t] / double(loc_u_new_n);
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

								for (l = 0; l < m_class_amount; ++l) {
									if (t == ofset * u_new_n)
										//new_weights[l] = mix_prob[l];
										new_weights[ofset * m_class_amount + l] = 1.0 / double(m_class_amount);
									m_prob_img[m_init_layer_idx[k] +
										(x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
										(y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
										l] = new_g_ij[t][l];
									/*if (m_layer_size[2 * k] == 4)
									{

										cout << m_init_layer_idx[k] << " " << k << endl;
										cout << x_min << " " << (t - ofset * u_new_n) << " " << y_l << endl;
										cout << y_min << " " << (y_min + (t - ofset * u_new_n)) << endl;
										cout << new_g_ij[t][l] << endl;
									}*/
									/*m_prob_img[m_init_layer_idx[k] +
									(x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
									(y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
									l] = 1.0/double(m_class_amount);*/
									//if (layer_size[k] == 32)
									//cout << t << " " << l << " " << new_g_ij[t][l] << endl;
									new_g_ij_0[t][l] = 0;
								}

							}
						}
					}
				}


			}
		}
//		else
//		{
////#pragma omp parallel
//			
//				int x_l = window_size, y_l = window_size, ofset = 0;
//				int itr, x_min, y_min, j, loc_u_new_n, t, l;
//				double pix_buf, cur_max, summ = 0, last_cur_max = 0, buf_max = 0;
//				bool stop_flag = true;
//				unsigned idx_max = 0;
//				for (l = ofset * m_class_amount; l < ofset * m_class_amount + m_class_amount; ++l)
//					new_weights[l] = 1.0 / double(m_class_amount);
////#pragma omp for
//				for (int r = 0; r < amount_window_x; ++r) {
//					x_min = r * window_size;
//					if (r < amount_window_x - 1)
//						x_l = window_size;
//					else
//						x_l = window_size + add_amount_x;
//
//					for (j = 0; j < amount_window_y; ++j) {
//						y_min = j * window_size;
//
//						if (j < amount_window_y - 1)
//							y_l = window_size;
//						else
//							y_l = window_size + add_amount_y;
//
//						itr = 0;
//						stop_flag = true;
//						cur_max = 0;
//						loc_u_new_n = y_l * x_l;
//
//						while (stop_flag && (itr < 500)) {
//							++itr;
//
//							for (l = ofset * u_new_n; l < ofset * u_new_n + loc_u_new_n; ++l) {
//								summ = 0;
//
//								pix_buf = m_image->get_image()[m_layer_idx[k] +
//									(x_min + (l - ofset * u_new_n) / y_l)*m_layer_size[2 * k + 1] + y_min + (l - ofset * u_new_n) % y_l];
//								for (t = 0; t < m_class_amount; ++t) {
//									if (m_image->get_mixture_type() == NORMAL)
//										summ += new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi))
//										*exp(-((pix_buf - mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
//									else {
//										if (m_image->get_mixture_type() == RAYLEIGH)
//											summ += new_weights[ofset * m_class_amount + t] * (pix_buf / (mix_scale[t] * mix_scale[t]))
//											* exp(-((pix_buf)*(pix_buf)) / (2.0 * mix_scale[t] * mix_scale[t]));
//
//										else {
//											if (m_image->get_mixture_type() == LOGNORMAL)
//												summ += new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*pix_buf))*exp(-((log(pix_buf)
//													- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
//										}
//									}
//									/*cout << "l  " << summ << " " << pix_buf << " " << layer_idx[k] +
//										(x_min + (l - ofset * u_new_n) / y_l)*layer_size[k] + y_min + (l - ofset * u_new_n) % y_l
//										<< " " << layer_idx[layer_amount-1]+1024*1024 << endl;*/
//
//								}
//
//								for (t = 0; t < m_class_amount; ++t) {
//									if (l == ofset * u_new_n)
//										buf_new_weights[ofset * m_class_amount + t] = 0;
//									if (m_image->get_mixture_type() == NORMAL)
//										new_g_ij[l][t] = new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*summ))
//										* exp(-((pix_buf - mix_shift[t]) * (pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
//
//									else {
//										if (m_image->get_mixture_type() == RAYLEIGH)
//											new_g_ij[l][t] = new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))
//											* exp(-(pix_buf * pix_buf) / (2.0 * mix_scale[t] * mix_scale[t]));
//
//										else {
//											if (m_image->get_mixture_type() == LOGNORMAL) {
//												//cout << "l  " << l << endl;
//												new_g_ij[l][t] = new_weights[ofset * m_class_amount + t] * (1 / (mix_scale[t] * sq_pi*summ*pix_buf))
//													* exp(-((log(pix_buf) - mix_shift[t]) * (log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
//
//											}
//										}
//									}
//									buf_new_weights[ofset * m_class_amount + t] += new_g_ij[l][t];
//									if (l == ofset * u_new_n + loc_u_new_n - 1)
//										new_weights[ofset * m_class_amount + t] = buf_new_weights[t] / double(loc_u_new_n);
//									if (cur_max < abs(new_g_ij[l][t] - new_g_ij_0[l][t]))
//										cur_max = abs(new_g_ij[l][t] - new_g_ij_0[l][t]);
//									new_g_ij_0[l][t] = new_g_ij[l][t];
//								}
//							}
//
//							if (stop_flag) {
//								if (cur_max != 0)
//									last_cur_max = cur_max;
//								(cur_max < accuracy) ? stop_flag = false : cur_max = 0;
//							}
//						}
////#pragma omp critical
//						
//							//if (layer_size[k] == 32)
//						cout << "jjjjj" << endl;
//
//							for (t = ofset * u_new_n; t < ofset * u_new_n + loc_u_new_n; ++t) {
//
//								for (l = 0; l < m_class_amount; ++l) {
//									if (t == ofset * u_new_n)
//										//new_weights[l] = mix_prob[l];
//										new_weights[ofset * m_class_amount + l] = 1.0 / double(m_class_amount);
//									m_prob_img[m_init_layer_idx[k] +
//										(x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
//										(y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
//										l] = new_g_ij[t][l];
//									if (m_layer_size[2 * k] == 4)
//									{
//
//										cout << m_init_layer_idx[k] << " " << k << endl;
//										cout << x_min + (t - ofset * u_new_n) <<  endl;
//										cout << y_min + (y_min + (t - ofset * u_new_n)) << endl;
//										cout << new_g_ij[t][l] << endl;
//										cout << m_prob_img[m_init_layer_idx[k] +
//											(x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
//											(y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
//											l] << endl;
//									}
//									/*m_prob_img[m_init_layer_idx[k] +
//									(x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
//									(y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
//									l] = 1.0/double(m_class_amount);*/
//									//if (layer_size[k] == 32)
//									//cout << t << " " << l << " " << new_g_ij[t][l] << endl;
//									new_g_ij_0[t][l] = 0;
//								}
//
//							}
//						
//					}
//				}
//
//
//			
//		}
        for (int t = 0; t < u_new_n * thr_nmb; ++t) {
            delete[] new_g_ij[t];
            delete[] new_g_ij_0[t];
        }

        delete[] new_g_ij;
        delete[] new_g_ij_0;
		delete[] new_weights;
		delete[] buf_new_weights;
        auto end1 = std::chrono::steady_clock::now();
        auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
        cout << "elapsed_ms1  " << elapsed_ms1.count() << endl;

    }
}


void initial_prob_img::generate_init_probs_max_apost_opMP_V2() {
    int window_size = m_init_window_size;
    float accuracy = 0.001;
    for (int i = 0; i < m_class_amount; ++i)
        mix_scale[i] = mix_scale[i] / (pow(2, m_layer_amount));

    for (int k = 0; k < m_layer_amount; ++k) {
        window_size = m_init_window_size;
        int add_amount_x = m_layer_size[2*k] % window_size;
        int add_amount_y = m_layer_size[2 * k + 1] % window_size;
        int amount_window_x = m_layer_size[2*k] / window_size;
        int amount_window_y = m_layer_size[2 * k + 1] / window_size;
        if (amount_window_x == 0) {
            amount_window_x = 1;
            amount_window_y = 1;
            add_amount_x = 0;
            add_amount_y = 0;
            window_size = m_layer_size[2 * k + 1];
        }
        int u_new_n = (window_size + add_amount_x) *(window_size + add_amount_y);
        int thr_nmb = 4;
        double* buf_mass = new double[u_new_n *thr_nmb];

        for (int i = 0; i < m_class_amount; ++i)
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
            int loc_hyp_cl_amount = m_class_amount;
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
                        buf_mass[l] = m_image->get_image()[m_layer_idx[k] +
                                (x_min + (l - ofset * u_new_n) / y_l)*m_layer_size[2 * k + 1] + y_min + (l - ofset * u_new_n) % y_l];
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

                                m_prob_img[m_init_layer_idx[k] +
                                        (x_min + (t - ofset * u_new_n) / y_l) * m_layer_size[2 * k + 1] * m_class_amount +
                                        (y_min + (t - ofset * u_new_n) % y_l) * m_class_amount +
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
