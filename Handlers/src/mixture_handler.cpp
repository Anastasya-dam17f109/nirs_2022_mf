
#include "mixture_handler.h"


mixture_handler::mixture_handler(shared_ptr < mix_img_obj> img, unsigned h_classes, double acc)
{
	accuracy      = acc;
	img_mask_list = img->get_mask_list();
	hyp_cl_amount= img->get_class_amount();
	raw_image    = img->get_raw_image();
	mixture_type = img->get_mixture_type();
	img_l_x      = img->get_image_len().first;
	img_l_y      = img->get_image_len().second;
	min_trg_size = img->get_min_targ_size();
	mix_shift    = img->get_shift();
	mix_scale    = img->get_scale();
	gen_mix_filename = img->get_filename();
	mix_weight   = shared_ptr<double[]>(new double[hyp_cl_amount]);
	window_size  = 16;
	cout << " mix partition by " << h_classes << " components:" << endl;

	mixture_inicalization();
	
	////detect_results();
	////BIC();
	
}

mixture_handler::mixture_handler(shared_ptr < mix_img_obj> img, unsigned h_classes, double acc, bool flag)
{
	accuracy = acc;
	img_mask_list = img->get_mask_list();
	hyp_cl_amount = img->get_class_amount();
	//raw_image = img->get_raw_image();
	//mixture_type = img->get_mixture_type();
	img_l_x = img->get_image_len().first;
	img_l_y = img->get_image_len().second;
	//min_trg_size = img->get_min_targ_size();
	//mix_shift = img->get_shift();
	//mix_scale = img->get_scale();
	gen_mix_filename = img->get_filename();
	//mix_weight = shared_ptr<double[]>(new double[hyp_cl_amount]);
	//window_size = 16;
	//cout << " mix partition by " << h_classes << " components:" << endl;

	mixture_inicalization();

}

void mixture_handler::get_classification_from_file(string filename)
{
	unsigned i, j;
	ifstream load_classification;
	load_classification.open(filename);
	//load_classification >> files_amount;
	// создаем массив-результат
	//class_flag = new unsigned*[img_l_x];
	for (i = 0; i < img_l_x; i++) {
		//class_flag[i] = new unsigned[img_l_y];
		for (j = 0; j < img_l_y; j++)
			load_classification >> class_flag[i][j] ;
	}
	load_classification.close();

}

mixture_handler::mixture_handler(shared_ptr <initial_prob_img> img, unsigned h_classes, double acc) {
	accuracy = acc;
	m_image = img;
	img_mask_list = img->get_m_image()->get_mask_list();
	hyp_cl_amount = img->get_class_amount();
	raw_image = img->get_m_image()->get_raw_image();
	mixture_type = img->get_m_image()->get_mixture_type();
	img_l_x = img->get_image_len().first;
	img_l_y = img->get_image_len().second;
	
	mix_shift = img->get_m_image()->get_shift();
	mix_scale = img->get_m_image()->get_scale();
	gen_mix_filename = img->get_m_image()->get_filename();
	mix_weight = shared_ptr<double[]>(new double[hyp_cl_amount]);
	window_size = 16;
	cout << " mix partition by " << h_classes << " components:" << endl;

	mixture_inicalization();

}

//выделение памяти под изображение-результат

void mixture_handler::mixture_inicalization() {
	unsigned i, j;
	// создаем массив-результат
	class_flag = new unsigned*[img_l_x];
	for (i = 0; i < img_l_x; i++) {
		class_flag[i] = new unsigned[img_l_y];
		for (j = 0; j < img_l_y; j++)
			class_flag[i][j] = 0;
	}
}

////
//
//void mixture_handler::SEMalgorithm_opMP_rayleigh_spec() {
//	int h_w = (window_size / 2);
//	int amount_window = img_l_x / h_w - 1;
//	re_comp_shifts = new long double*[amount_window*amount_window];
//	re_comp_scales = new long double*[amount_window*amount_window];
//	bool     stop_flag = true;
//
//	int i_new_n = window_size * window_size;
//	int i_n_add = img_l_x - (img_l_x / h_w)* h_w;
//	int est_amount = 2;
//	unsigned i, t, r;
//	double        cur_max = 0;
//	const double  sq_pi = sqrt(2 * pi);
//	double        min_p = double(min_trg_size* min_trg_size) / double(window_size * window_size);
//
//	int x_min, y_min, i_lim, i_step, itr, horiz, vert, j;
//	int u_new_n = 0;
//	u_new_n = window_size * window_size;
//	horiz = window_size;
//	vert = window_size;
//	double d_u_new_n = double(u_new_n);
//
//	for (i = 0; i < amount_window*amount_window; ++i) {
//		re_comp_shifts[i] = new long double[hyp_cl_amount];
//		re_comp_scales[i] = new long double[hyp_cl_amount];
//		for (t = 0; t < hyp_cl_amount; ++t) {
//			re_comp_shifts[i][t] = 0;
//			re_comp_scales[i][t] = 0;
//		}
//	}
//
//
//	double   max_buf = 0;
//	double bic_value = 0;
//	double pre_bic_value = 0;
//	int pre_amount = 1;
//	auto begin1 = std::chrono::steady_clock::now();
//
//	double** new_g_ij = new double*[(window_size + i_n_add)*(window_size + i_n_add)];
//	double** new_g_ij_0 = new double*[(window_size + i_n_add)*(window_size + i_n_add)];
//	double**    y_i_j = new double*[(window_size + i_n_add)*(window_size + i_n_add)];
//
//	double * new_weights_mean = new  double[hyp_cl_amount];
//
//	double ** new_shifts = new double*[est_amount];
//	double ** new_scales = new double*[est_amount];
//
//	for (j = 0; j < est_amount; ++j) {
//		new_shifts[j] = new  double[hyp_cl_amount];
//		new_scales[j] = new  double[hyp_cl_amount];
//	}
//
//
//	double** new_med_sample = new double*[hyp_cl_amount];
//
//	double * pre_new_weights = new double[hyp_cl_amount];
//	double * pre_new_shifts = new double[hyp_cl_amount];
//	double * pre_new_scales = new double[hyp_cl_amount];
//
//	double * max_L_mass = new  double[est_amount];
//
//	for (j = 0; j < hyp_cl_amount; ++j) {
//		pre_new_weights[j] = 0;
//		pre_new_shifts[j] = 0;
//		pre_new_scales[j] = 0;
//		new_med_sample[j] = new double[(window_size + i_n_add)*(window_size + i_n_add)];
//
//		new_weights_mean[j] = 0;
//
//		for (int l = 0; l < est_amount; l++) {
//			new_shifts[l][j] = 0;
//			new_scales[l][j] = 0;
//			max_L_mass[l] = 0;
//		}
//	}
//
//	double val1, psumm1;
//	boost::random::uniform_01 <> dist_poly1;
//	boost::random::mt19937 generator1{ static_cast<std::uint32_t>(time(0)) };
//	for (j = 0; j < (window_size + i_n_add)*(window_size + i_n_add); ++j) {
//		new_g_ij[j] = new double[hyp_cl_amount];
//		new_g_ij_0[j] = new double[hyp_cl_amount];
//		y_i_j[j] = new double[hyp_cl_amount];
//	}
//	bool* med_flag = new bool[est_amount];
//	med_flag[0] = true;
//	for (int r = 0; r < amount_window; ++r) {
//		x_min = r * h_w;
//		for (j = 0; j < amount_window; ++j) {
//			y_min = j * h_w;
//			for (int itr_cl_am = 1; itr_cl_am <= hyp_cl_amount; ++itr_cl_am) {
//				stop_flag = true;
//				itr = 0;
//				/*for (int l = 0; l < (window_size + i_n_add)*(window_size + i_n_add); ++l) {
//					psumm1 = 0;
//					for (int t = 0; t < itr_cl_am; t++) {
//						val1 = dist_poly1(generator1);
//						psumm1 += val1;
//						new_g_ij[l][t] = val1;
//						new_g_ij_0[l][t] = 0;
//					}
//					for (int t = 0; t < itr_cl_am; t++)
//						new_g_ij[l][t] = new_g_ij[l][t] / psumm1;
//				}
//*/
//#pragma omp parallel 
//				{
//					double val, psumm;
//					int l;
//					int u_new_n = window_size * window_size;
//					int loc_hyp_cl_amount = itr_cl_am;
//					boost::random::uniform_01 <> dist_poly;
//					boost::random::mt19937 generator{ static_cast<std::uint32_t>(time(0)) };
//#pragma omp for
//					for (int t = 0; t < u_new_n; ++t) {
//						psumm = 0;
//						for (l = 0; l < loc_hyp_cl_amount; ++l) {
//							val = dist_poly(generator);
//							psumm += val;
//							new_g_ij[t][l] = val;
//							new_g_ij_0[t][l] = 0;
//						}
//
//						for (l = 0; l < loc_hyp_cl_amount; ++l)
//							new_g_ij[t][l] = new_g_ij[t][l] / psumm;
//					}
//				}
//				/*
//					if (r != amount_window - 1)
//						horiz = window_size;
//					else
//						horiz = window_size + i_n_add;
//					if (j != amount_window - 1)
//						vert = window_size;
//					else
//						vert = (window_size + i_n_add);
//					u_new_n = horiz * vert;
//				*/
//
//				while (stop_flag &&  itr < 300) {
//					++itr;
//					cur_max = 0;
//#pragma omp parallel 
//					{
//						int loc_hyp_cl_amount = itr_cl_am;
//						string mixture_type = mixture_type;
//						int u_new_n = window_size * window_size;
//						int t;
//						int loc_horiz = horiz;
//						int loc_x_min = x_min;
//						int loc_y_min = y_min;
//						double bound_d, bound_u, psumm, val, pix_buf;
//						boost::random::uniform_01 <> dist_poly;
//						boost::random::mt19937 generator{ static_cast<std::uint32_t>(time(0)) };
//#pragma omp for  
//						for (int l = 0; l < u_new_n; ++l) {
//							pix_buf = my_picture[loc_x_min + l / loc_horiz][loc_y_min + l % loc_horiz];
//							bound_d = 0;
//							bound_u = new_g_ij[l][0];
//							val = dist_poly(generator);
//							for (t = 0; t < loc_hyp_cl_amount; ++t) {
//								if ((val < bound_u)
//									&& (val >= bound_d)) {
//									y_i_j[l][t] = 1;
//# pragma omp critical
//									{
//										new_med_sample[t][int(new_weights_mean[t])] = pix_buf;
//										++new_weights_mean[t];
//										new_scales[0][t] += pix_buf * pix_buf;
//										
//									}
//
//								}
//								else
//									y_i_j[l][t] = 0;
//								bound_d += new_g_ij[l][t];
//								if (t < loc_hyp_cl_amount - 2)
//									bound_u += new_g_ij[l][t + 1];
//								else
//									if (t == loc_hyp_cl_amount - 2)
//										bound_u = 1;
//							}
//						}
//					}
//
//
//					med_flag[1] = true;
//					for (int l = 0; l < itr_cl_am; ++l) {
//							new_shifts[0][l] = 0;
//							new_shifts[1][l] = 0;
//							new_scales[0][l] = sqrt(new_scales[0][l] / double(2 * new_weights_mean[l]));
//
//							new_scales[1][l] = find_med(new_med_sample[l], new_weights_mean[l]) / sqrt(log(4.0));
//							if (new_scales[1][l] < 0)
//								med_flag[1] = false;
//					}
//#pragma omp parallel 
//					{
//						
//						double   pix_buf = 0; int l;
//						int loc_horiz = horiz;
//						int loc_x_min = x_min;
//						int loc_y_min = y_min;
//						int loc_hyp_cl_amount = itr_cl_am;
//						int u_new_n = window_size * window_size;
//						double d_u_new_n = double(u_new_n);
//
//#pragma omp for  
//						for (int t = 0; t < u_new_n; ++t) {
//							pix_buf = my_picture[loc_x_min + t / loc_horiz][loc_y_min + t % loc_horiz];
//							for (l = 0; l < loc_hyp_cl_amount; ++l) {
//
//# pragma omp critical
//								{
//									
//									if (cur_max < abs(new_g_ij[t][l] - new_g_ij_0[t][l]))
//										cur_max = abs(new_g_ij[t][l] - new_g_ij_0[t][l]);
//								}
//								new_g_ij_0[t][l] = new_g_ij[t][l];
//							}
//						}
//					}
//
//
//					for (int l = 0; l < itr_cl_am; ++l) {
//						
//						new_weights_mean[l] = new_weights_mean[l] / d_u_new_n;
//					}
//#pragma omp parallel
//					{
//
//						double B, pix_buf;
//						
//						bool pre_flag = true;
//						bool flag = false;
//						int l, t;
//						int u_new_n = window_size * window_size;
//						int loc_x_min = x_min;
//						int loc_y_min = y_min;
//						int loc_est_amount = est_amount;
//						int loc_horiz = horiz;
//						int loc_itr_cl_am = itr_cl_am;
//						double d_u_new_n = double(u_new_n);
//#pragma omp for
//						for (int m = 0; m < u_new_n; m++) {
//							pix_buf = my_picture[loc_x_min + m / loc_horiz][loc_y_min + m % loc_horiz];
//							flag = false;
//							for (l = 0; l < loc_est_amount; ++l) {
//								B = 0;
//								for (t = 0; t < loc_itr_cl_am; ++t) {
//									if (((new_scales[l][t] != 0) && (new_weights_mean[t] != 0)) && (med_flag[l])) {
//											B += new_weights_mean[t] * (pix_buf / (new_scales[l][t] * new_scales[l][t]))*
//											exp(-(pow(pix_buf, 2)) /
//											(2.0 * new_scales[l][t] * new_scales[l][t]));
//									}
//									else {
//										flag = true;
//										break;
//									}
//								}
//#pragma omp critical
//								if (B > 0)
//									max_L_mass[l] += log(B / (d_u_new_n));
//								else
//									max_L_mass[l] += 0;
//							}
//						}
//
//					}
//
//					int max_idx = 0;
//
//					for (int m = 0; m < est_amount; ++m) {
//						//cout << m << " " << max_L_mass[m] << endl;
//						if ((max_L_mass[max_idx] < max_L_mass[m]) && ((med_flag[m])))
//							max_idx = m;
//					}
//
//					for (int m = 0; m < itr_cl_am; ++m) {
//						new_shifts[0][m] = new_shifts[max_idx][m];
//						new_scales[0][m] = new_scales[max_idx][m];
//					}
//					//cout << "max_idx " << max_idx<<" "<< new_shifts[1][0]<<" "<< new_scales[1][0]<< " " << new_shifts[1][1] << " " << new_scales[1][1] << " "<< med_flag[1]<<endl;
//					//cout << "max_idx_t "  << " " << new_shifts[0][0] << " " << new_scales[0][0] << " " << new_shifts[0][1] << " " << new_scales[0][1] << " " << med_flag[1] << endl;
//#pragma omp parallel 
//					{
//						double   pix_buf = 0;
//						double summ;
//						int l;
//					
//						int loc_horiz = horiz;
//						int loc_x_min = x_min;
//						int loc_y_min = y_min;
//						int loc_hyp_cl_amount = itr_cl_am;
//						int u_new_n = window_size * window_size;
//						double d_u_new_n = double(u_new_n);
//#pragma omp for  
//						for (int t = 0; t < u_new_n; ++t) {
//							pix_buf = my_picture[loc_x_min + t / loc_horiz][loc_y_min + t % loc_horiz];
//							summ = 0;
//							for (l = 0; l < loc_hyp_cl_amount; ++l) {
//								if ((new_weights_mean[l] >= 0)
//									&& (new_scales[0][l] > 0)) {
//									
//											summ += new_weights_mean[l] * (pix_buf / (new_scales[0][l] * new_scales[0][l]))*exp(-((pix_buf
//												)*(pix_buf)) / (2.0 * new_scales[0][l] * new_scales[0][l]));
//									
//								}
//
//							}
//							for (l = 0; l < loc_hyp_cl_amount; ++l)
//								if ((new_weights_mean[l] >= 0)
//									&& (new_scales[0][l] > 0))
//									
//											new_g_ij[t][l] = new_weights_mean[l] * (pix_buf / (new_scales[0][l] * new_scales[0][l]* summ))*exp(-(pix_buf
//												*pix_buf) / (2.0 * new_scales[0][l] * new_scales[0][l]));
//									
//								else
//									new_g_ij[t][l] = 0;
//
//						}
//					}
//
//					if (stop_flag) {
//						if (cur_max < accuracy)
//							stop_flag = false;
//						else {
//							cur_max = 0;
//							for (int l = 0; l < itr_cl_am; ++l) {
//								new_weights_mean[l] = 0;
//								for (int t = 0; t < est_amount; ++t) {
//									new_shifts[t][l] = 0;
//									new_scales[t][l] = 0;
//									max_L_mass[t] = 0;
//								}
//							}
//						}
//					}
//
//				}
//
//				for (int l = 0; l < itr_cl_am; ++l) {
//					new_weights_mean[l] = 0;
//					//cout<<"dd "<< new_shifts[0][l]<< " "<<new_scales[0][l] <<endl;
//				}
//#pragma omp parallel 
//				{
//					double buf_max = 0;
//					int      idx_max = 0;
//					double val, psumm;
//					int l;
//					int loc_hyp_cl_amount = itr_cl_am;
//					boost::random::uniform_01 <> dist_poly;
//					int u_new_n = window_size * window_size;
//					boost::random::mt19937 generator{ static_cast<std::uint32_t>(time(0)) };
//#pragma omp for
//					for (int t = 0; t < u_new_n; ++t) {
//
//						buf_max = new_g_ij[t][0];
//						idx_max = 0;
//						psumm = 0;
//						for (l = 0; l < loc_hyp_cl_amount; ++l) {
//							if (buf_max < new_g_ij[t][l]) {
//								buf_max = new_g_ij[t][l];
//								idx_max = l;
//							}
//							new_g_ij_0[t][l] = 0;
//							val = dist_poly(generator);
//							psumm += val;
//							new_g_ij[t][l] = val;
//							new_g_ij_0[t][l] = 0;
//						}
//#pragma omp critical
//						++new_weights_mean[idx_max];
//						for (l = 0; l < loc_hyp_cl_amount; ++l)
//							new_g_ij[t][l] = new_g_ij[t][l] / psumm;
//						//new_g_ij[t][l] = 1.0 / double(k);
//
//					}
//				}
//				double big_summ = 0;
//#pragma omp parallel
//				{
//					string mixture_type = mixture_type;
//					double pix_buf, summ;
//					int t; int loc_horiz = horiz;
//					int loc_x_min = x_min;
//					int loc_y_min = y_min;
//#pragma omp for
//					for (int l = 0; l < u_new_n; ++l) {
//						summ = 0;
//						pix_buf = my_picture[loc_x_min + l / loc_horiz][loc_y_min + l % loc_horiz];
//						for (t = 0; t < itr_cl_am; ++t)
//							
//									summ += new_weights_mean[t] * (pix_buf / (new_scales[0][t] * new_scales[0][t]))*exp(-((pix_buf
//										)*(pix_buf)) / (2.0 * new_scales[0][t] * new_scales[0][t]));
//							
//# pragma omp critical
//						big_summ += log(summ);
//					}
//
//				}
//				bic_value = 2 * big_summ - log(u_new_n)*(2 * itr_cl_am - 1);
//
//				if (itr_cl_am == 1) {
//					pre_bic_value = bic_value;
//					pre_amount = 1;
//					for (int l = 0; l < itr_cl_am; ++l) {
//						pre_new_weights[l] = new_weights_mean[l];
//						pre_new_shifts[l] = new_shifts[0][l];
//						pre_new_scales[l] = new_scales[0][l];
//					}
//				}
//				else {
//					if (bic_value > pre_bic_value) {
//						pre_bic_value = bic_value;
//						pre_amount = itr_cl_am;
//						for (int l = 0; l < itr_cl_am; ++l) {
//							pre_new_weights[l] = new_weights_mean[l];
//							pre_new_shifts[l] = new_shifts[0][l];
//							pre_new_scales[l] = new_scales[0][l];
//						}
//					}
//				}
//				for (int l = 0; l < itr_cl_am; ++l) {
//
//					new_weights_mean[l] = 0;
//					for (int t = 0; t < est_amount; ++t) {
//						new_shifts[t][l] = 0;
//						new_scales[t][l] = 0;
//					}
//				}
//
//			}
//			//cout << r << " " << j << " " << pre_amount<<endl;
//			for (int l = 0; l < pre_amount; ++l) {
//				//if(pre_amount==2)
//				//cout <<"pre_new_scales[l] "<< pre_new_scales[l] << endl;
//				if ((pre_new_weights[l] / d_u_new_n > min_p)
//					&& (pre_new_scales[l] > 0)) {
//
//					re_comp_shifts[r *  amount_window + j][l] = pre_new_scales[l];
//					re_comp_scales[r * amount_window + j][l] = pre_new_scales[l];
//				}
//				pre_new_weights[l] = 0;
//				pre_new_shifts[l] = 0;
//				pre_new_scales[l] = 0;
//			}
//		}
//	}
//
//	auto end1 = std::chrono::steady_clock::now();
//	auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
//	cout << "elapsed_ms1  " << elapsed_ms1.count() << "\n";
//	for (int t = 0; t < (window_size + i_n_add)*(window_size + i_n_add); ++t) {
//		delete[] new_g_ij[t];
//		delete[] y_i_j[t];
//		delete[] new_g_ij_0[t];
//	}
//
//	delete[] new_g_ij;
//	delete[] y_i_j;
//	delete[] new_g_ij_0;
//	delete[] med_flag;
//	delete[] new_weights_mean;
//	delete[] pre_new_weights;
//	delete[] pre_new_shifts;
//	delete[] pre_new_scales;
//	delete[] max_L_mass;
//
//	for (j = 0; j < est_amount; ++j) {
//		delete[] new_shifts[j];
//		delete[] new_scales[j];
//	}
//	delete[] new_shifts;
//	delete[] new_scales;
//
//	for (j = 0; j < hyp_cl_amount; ++j) {
//		delete[] new_med_sample[j];
//
//	}
//	delete[] new_med_sample;
//
//}

//раскраска картики - em алгоритм , openMP version
// иначе - разделение картинки сеточным методом разделения смесей 

void mixture_handler::mixture_optimal_redraw_opMP() {

	int add_amount_x = img_l_x% window_size;
	int add_amount_y = img_l_y % window_size;
	int amount_window_x = img_l_x / window_size ;
	int amount_window_y = img_l_y / window_size;
	unsigned u_new_n = (window_size + add_amount_x) *(window_size + add_amount_y);

	
	auto begin1 = std::chrono::steady_clock::now();
	#pragma omp parallel
	{
		int loc_window_size = window_size;
		int loc_hyp_cl_amount = hyp_cl_amount;
		int x_l = loc_window_size;
		int y_l = loc_window_size;

		int itr, x_min, y_min, j, loc_u_new_n, t, l;
		double** new_g_ij = new double*[u_new_n];
		double** new_g_ij_0 = new double*[u_new_n];
		double * new_weights = new  double[loc_hyp_cl_amount];
		double * new_shifts = new  double[loc_hyp_cl_amount];
		double * new_scales = new  double[loc_hyp_cl_amount];
		double * buf_new_weights = new  double[loc_hyp_cl_amount];
		
		for (l = 0; l < u_new_n; ++l) {
			new_g_ij[l] = new double[loc_hyp_cl_amount];
			new_g_ij_0[l] = new double[loc_hyp_cl_amount];
			for (t = 0; t < loc_hyp_cl_amount; t++) {
				new_g_ij[l][t] = 0;
				new_g_ij_0[l][t] = 0;
			}
		}
		for (l = 0; l < loc_hyp_cl_amount; ++l) {
			new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
			new_shifts[l] = mix_shift[l];
			new_scales[l] = mix_scale[l];
		}
		
		double summ = 0;
		const double sq_pi = sqrt(2 * pi);
		double pix_buf,cur_max;
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
					
					for (l = 0; l < loc_u_new_n; ++l) {
						summ = 0;
						pix_buf = raw_image[(x_min + l / y_l)*img_l_x + y_min + l % y_l];
						for (t = 0; t < loc_hyp_cl_amount; ++t) {
						//if (mixture_type == "normal")
							//summ += sum_func(mix_shift[t], mix_scale[t], new_weights[t], pix_buf);
							//cout << sum_func(mix_shift[t], mix_scale[t], new_weights[t], pix_buf)<< endl;
								/*summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi))*exp(-((pix_buf
									- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));*/
							/*else {
								if (mixture_type == "rayleigh")*/
									/*summ += new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-((pix_buf
										)*(pix_buf )) / (2.0 * mix_scale[t] * mix_scale[t]));*/
							/*}*/
							//else {
							//if (mixture_type == "lognormal")
								summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi*pix_buf))*exp(-((log(pix_buf)
									- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
						//}

						}

						for (t = 0; t < loc_hyp_cl_amount; ++t) {
								if (l == 0)
									buf_new_weights[t] = 0;
								//if (mixture_type == "normal")
									/*new_g_ij[l][t] = new_weights[t] * (1 / (mix_scale[t] * sq_pi*summ))*exp(-((pix_buf
										- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));*/
								//new_g_ij[l][t] = sum_func(mix_shift[t], mix_scale[t], new_weights[t], pix_buf) / summ; 
								/*else {
									if (mixture_type == "rayleigh")*/
										/*new_g_ij[l][t] = new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-(pix_buf
											*pix_buf) / (2.0 * mix_scale[t] * mix_scale[t]));*/
								/*}*/
								///else {
								//if (mixture_type == "lognormal")
									new_g_ij[l][t] = new_weights[t] * (1 / (mix_scale[t] * sq_pi*summ*pix_buf))*exp(-((log(pix_buf)
										- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
								//}
								buf_new_weights[t] += new_g_ij[l][t];
								if (l == loc_u_new_n - 1)
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
					//cout << "jjjjj" << endl;

					for (t = 0; t < loc_u_new_n; ++t) {
						buf_max = new_g_ij[t][0];
						idx_max = 0;
						for (l = 0; l < loc_hyp_cl_amount; ++l) {
							if (t == 0)
								//new_weights[l] = mix_prob[l];
								new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
							if (buf_max < new_g_ij[t][l]) {
								buf_max = new_g_ij[t][l];
								idx_max = l;
							}
							new_g_ij_0[t][l] = 0;
						}
						class_flag[x_min + t / y_l][y_min + t % y_l] = idx_max + 1;
					}


				}
			}

		}

		for (t = 0; t < u_new_n; ++t) {
			delete[] new_g_ij[t];
			delete[] new_g_ij_0[t];
		}
		delete[] new_g_ij;
		delete[] new_g_ij_0;
		delete[] new_weights;
		delete[] buf_new_weights;
		delete[] new_shifts;
		delete[] new_scales;
	}
	
	auto end1 = std::chrono::steady_clock::now();
	auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
	cout << "elapsed_ms1  " << elapsed_ms1.count() <<endl;
}

// попытка оптимизировать верхнюю функцию

void mixture_handler::mixture_optimal_redraw_opMP_V2(){

	int add_amount_x = img_l_x % window_size;
	int add_amount_y = img_l_y % window_size;
	int amount_window_x = img_l_x / window_size;
	int amount_window_y = img_l_y / window_size;
	int u_new_n = (window_size + add_amount_x) *(window_size + add_amount_y);
	int thr_nmb = 4;
	double** new_g_ij = new double*[u_new_n *thr_nmb];
	double** new_g_ij_0 = new double*[u_new_n * thr_nmb];
	cout << omp_get_num_threads() << endl;
	for (int l = 0; l < u_new_n * thr_nmb; ++l) {
		new_g_ij[l] = new double[hyp_cl_amount];
		new_g_ij_0[l] = new double[hyp_cl_amount];
		for (int t = 0; t < hyp_cl_amount; t++) {
			new_g_ij[l][t] = 0;
			new_g_ij_0[l][t] = 0;
		}
	}
	auto begin1 = std::chrono::steady_clock::now();
#pragma omp parallel
	{
		int loc_window_size = window_size;
		int loc_hyp_cl_amount = hyp_cl_amount;
		int x_l = loc_window_size;
		int y_l = loc_window_size;
		int ofset = omp_get_thread_num();
		//cout << "ofset  " << ofset << endl;
		int itr, x_min, y_min, j, loc_u_new_n, t, l;

		double * new_weights = new  double[loc_hyp_cl_amount];
		double * new_shifts = new  double[loc_hyp_cl_amount];
		double * new_scales = new  double[loc_hyp_cl_amount];
		double * buf_new_weights = new  double[loc_hyp_cl_amount];


		for (l = 0; l < loc_hyp_cl_amount; ++l) {
			new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
			new_shifts[l] = mix_shift[l];
			new_scales[l] = mix_scale[l];
		}

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
						pix_buf = raw_image[(x_min + (l - ofset * u_new_n) / y_l)*img_l_x + y_min + (l - ofset * u_new_n) % y_l];
						for (t = 0; t < loc_hyp_cl_amount; ++t) {
							if (mixture_type == NORMAL)
								summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi))*exp(-((pix_buf
									- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
							else {
								if (mixture_type == RAYLEIGH)
									summ += new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-((pix_buf
										)*(pix_buf)) / (2.0 * mix_scale[t] * mix_scale[t]));

								else {
									if (mixture_type == LOGNORMAL)
							summ += new_weights[t] * (1 / (mix_scale[t] * sq_pi*pix_buf))*exp(-((log(pix_buf)
								- mix_shift[t])*(log(pix_buf) - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));
								}
							}

						}

						for (t = 0; t < loc_hyp_cl_amount; ++t) {
							if (l == ofset * u_new_n)
								buf_new_weights[t] = 0;
							if (mixture_type == NORMAL)
								new_g_ij[l][t] = new_weights[t] * (1 / (mix_scale[t] * sq_pi*summ))*exp(-((pix_buf
									- mix_shift[t])*(pix_buf - mix_shift[t])) / (2.0 * mix_scale[t] * mix_scale[t]));

							else {
								if (mixture_type == RAYLEIGH)
									new_g_ij[l][t] = new_weights[t] * (pix_buf / (mix_scale[t] * mix_scale[t]))*exp(-(pix_buf
										*pix_buf) / (2.0 * mix_scale[t] * mix_scale[t]));

								else {
									if (mixture_type == LOGNORMAL) {
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
					//cout << "jjjjj" << endl;

					for (t = ofset * u_new_n; t < ofset * u_new_n + loc_u_new_n; ++t) {
						buf_max = new_g_ij[t][0];
						idx_max = 0;
						for (l = 0; l < loc_hyp_cl_amount; ++l) {
							if (t == ofset * u_new_n)
								//new_weights[l] = mix_prob[l];
								new_weights[l] = 1.0 / double(loc_hyp_cl_amount);
							if (buf_max < new_g_ij[t][l]) {
								buf_max = new_g_ij[t][l];
								idx_max = l;
							}
							new_g_ij_0[t][l] = 0;
						}
						class_flag[x_min + (t - ofset * u_new_n) / y_l][y_min + (t - ofset * u_new_n) % y_l] = idx_max + 1;
					}
				}
			}
		}
		delete[] new_weights;
		delete[] buf_new_weights;
		delete[] new_shifts;
		delete[] new_scales;
	}
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

//раскраска картики - по локальным областям, с использованием квадродеревьев
//, openMP version

void mixture_handler::q_tree_optimal_redraw_opMP(){

    int amount_window_x = img_l_x / window_size ;
    int amount_window_y = img_l_y / window_size;


    auto begin1 = std::chrono::steady_clock::now();
    #pragma omp parallel
    {
       
        quad_tree_handler tree = quad_tree_handler(m_image, window_size, class_flag,4);

        #pragma omp for
        for (int i = 0; i < amount_window_x*amount_window_y; ++i) {
            
                tree.set_probabilities(i/ amount_window_x, i%amount_window_x);
                tree.bottom_up_pass();
				tree.up_down_pass_V2();
            // tree.up_down_pass();
                // возможно, на сюда придется секцию critical
#pragma omp critical
				{
					tree.split_image_by_vote();
				}
            
        }
  }

    auto end1 = std::chrono::steady_clock::now();
    auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
    cout << "elapsed_ms1  " << elapsed_ms1.count() <<endl;
}

// раскраска по колмогорову, с распараллеливанием

void mixture_handler::kolmogorov_optimal_redraw_opMP() {
	int length_ = 10;
	auto begin1 = std::chrono::steady_clock::now();
#pragma omp parallel
	{
		int half_step = length_ / 2, j , k, l, m;
		int loc_hyp_cl_amount = hyp_cl_amount;
		int iters_i_x = img_l_x / half_step;
		int iters_i_y = img_l_y / half_step;
		while ((iters_i_x*half_step + length_) > img_l_x)
			iters_i_x--;
		while ((iters_i_y*half_step + length_) > img_l_y)
			iters_i_y--;
		int x_l, y_l, idx_class;
		double* buf_img = new double[4 * length_*length_];

		// вычисление максимума правдоподобия
		auto L_max_calculation = [&](double* data, int data_size, int iter, double mix_shift, double mix_scale, double* max_L_mass) {
			double buf_max_l = 0;
			bool flag = false;
			double B;
			for (int m = 0; m < data_size; m++) {
				B = 0;
				if (mix_scale != 0) {
					if (mixture_type == NORMAL)
						B = (1.0 / (mix_scale))*
							exp(-(pow(data[m] - mix_shift, 2)) /
							(2.0 * mix_scale * mix_scale));
					else 
						if (mixture_type == LOGNORMAL)
							B = (1.0 / (mix_scale*data[m]))*
								exp(-(pow(log(data[m]) - mix_shift, 2)) /
								(2.0 * mix_scale * mix_scale));
						else
							if (mixture_type == RAYLEIGH)
								B = (data[m] / pow(mix_scale, 2))*
									exp(-(pow(data[m], 2)) /
									(2.0 * mix_scale * mix_scale));
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

		// классификация по значениям статистики колмогорова по интервалам
		auto kolmogorov_stats = [&](double* data, int data_size, 
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			unsigned i, j;
			int intervals_amount = 60;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr;
			int* flag_mass = new int[hyp_cl_amount];
			double* stats_mass = new double[hyp_cl_amount];
			int flag_summ = 0;
			int flag_idx = 0;
			double* nu_i = new double[intervals_amount];
			double max_value = find_k_stat(data, data_size, data_size - 1) + 1;
			double min_value = find_k_stat(data, data_size, 0);
			double len_interval = (max_value - min_value) / intervals_amount;
			double *max_L_mass = new double[hyp_cl_amount];
			for (i = 0; i < intervals_amount; ++i)
				nu_i[i] = 0;
			for (i = 0; i < data_size; ++i) {
				k = 0;
				flag = true;
				while (flag) {
					if (((len_interval * k + min_value) <= data[i]) && ((len_interval * (k + 1) + min_value) > data[i]))
						flag = false;
					else
						k++;
				}
				nu_i[k] = nu_i[k] + 1;
			}

			double dn_bound;
			if (length_ == 5)
				dn_bound = 0.264;
			else
				dn_bound = 0.134;
			//double dn_bound = 0.238;
			for (i = 0; i < hyp_cl_amount; ++i) {
				F_n_curr = 0;
				max_d_n = 0;
				if (mixture_type == NORMAL)
					max_d_n = cdf(normal(mix_shift[i], mix_scale[i]), min_value);
				if (mixture_type == RAYLEIGH)
					max_d_n = cdf(rayleigh(mix_scale[i]), min_value);
				for (j = 0; j < intervals_amount; ++j) {
					F_n_curr += nu_i[j] / data_size;
					if (j != intervals_amount - 1) {
						for (k = 0; k < 2; ++k) {
							if (mixture_type == NORMAL)
								buf_d_n = abs(cdf(normal(mix_shift[i], mix_scale[i]), (len_interval *(j + k) + min_value)) - F_n_curr);
							if (mixture_type == RAYLEIGH)
								buf_d_n = abs(cdf(rayleigh(mix_scale[i]), (len_interval *(j + k) + min_value)) - F_n_curr);
							if (buf_d_n > max_d_n)
								max_d_n = buf_d_n;
						}
					}
					else {
						if (mixture_type == NORMAL)
							buf_d_n = abs(cdf(normal(mix_shift[i], mix_scale[i]), (len_interval *(j + 1) + min_value)) - 1);
						if (mixture_type == RAYLEIGH)
							buf_d_n = abs(cdf(rayleigh(mix_scale[i]), (len_interval *(j + 1) + min_value)) - 1);
						if (buf_d_n > max_d_n)
							max_d_n = buf_d_n;
					}
				}

				stats_mass[i] = max_d_n;
				if (max_d_n < dn_bound)
					flag_mass[i] = 1;
				else
					flag_mass[i] = 0;
			}
			for (i = 0; i < hyp_cl_amount; ++i) {
				flag_summ += flag_mass[i];
				if (flag_mass[i] == 1)
					flag_idx = i;
			}
			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] stats_mass;
			delete[] max_L_mass;
			delete[] nu_i;
			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};

		// классификация по значениям статистики максимального правдоподобия
		auto max_likehood_stats = [&](double* data, int data_size, 
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			unsigned i, j;
			int intervals_amount = 60;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr;
			int flag_summ = 0;
			int flag_idx = 0;
			double *max_L_mass = new double[hyp_cl_amount];
			for (int i = 0; i < hyp_cl_amount; ++i)
				L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

			int beg_idx = 0;
			int max_idx = 0;

			for (int m = beg_idx; m < hyp_cl_amount; ++m) 
				if (max_L_mass[max_idx] < max_L_mass[m])
					max_idx = m;
			
			flag_summ = 1;
			flag_idx = max_idx;


			delete[] max_L_mass;

			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};

		// классификация по значениям статистики колмогорова
		auto kolmogorov_stats2 = [&](double* data, int data_size, 
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			unsigned i, j;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr;
			int* flag_mass = new int[hyp_cl_amount];
			double* stats_mass = new double[hyp_cl_amount];
			int flag_summ = 0;
			int flag_idx = 0;
			double *max_L_mass = new double[hyp_cl_amount];
			quickSort(data, 0, data_size - 1, data_size / 2);
			double dn_bound;
			if (length_ == 5)
				dn_bound = 0.264;
			else
				dn_bound = 0.134;
			//double dn_bound = 0.238;
			for (i = 0; i < hyp_cl_amount; ++i) {
				F_n_curr = 0;
				max_d_n = 0;
				if (mixture_type == LOGNORMAL)
					max_d_n = cdf(lognormal(mix_shift[i], mix_scale[i]), data[0]);
				else {
					if (mixture_type == NORMAL)
						max_d_n = cdf(normal(mix_shift[i], mix_scale[i]), data[0]);
					else
						if (mixture_type == RAYLEIGH)
							max_d_n = cdf(rayleigh(mix_scale[i]), data[0]);
				}
				
				for (j = 1; j < data_size; ++j) {
					F_n_curr += 1.0 / data_size;
					if (j != data_size - 1) {
						for (k = 0; k < 2; ++k) {
							if (mixture_type == LOGNORMAL)
								buf_d_n = abs(cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]) - (F_n_curr - k * 1.0 / data_size));
							else {
								if (mixture_type == NORMAL)
									buf_d_n = abs(cdf(normal(mix_shift[i], mix_scale[i]), data[j]) - (F_n_curr - k * 1.0 / data_size));
								else {
									if (mixture_type == RAYLEIGH)
										buf_d_n = abs(cdf(rayleigh(mix_scale[i]), data[j]) - (F_n_curr - k * 1.0 / data_size));
								}
							}
						}
					}
					else {
						if (mixture_type == LOGNORMAL)
							buf_d_n = abs(cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]) - 1);
						else {
							if (mixture_type == NORMAL)
								buf_d_n = abs(cdf(normal(mix_shift[i], mix_scale[i]), data[j]) - 1);
							else {
								if (mixture_type == RAYLEIGH)
									buf_d_n = abs(cdf(rayleigh(mix_scale[i]), data[j]) - 1);
							}
						}	
					}
					if (buf_d_n > max_d_n)
						max_d_n = buf_d_n;
				}
				
				stats_mass[i] = max_d_n;
				if (max_d_n < dn_bound)
					flag_mass[i] = 1;
				else
					flag_mass[i] = 0;
			}
			for (i = 0; i < hyp_cl_amount; ++i) {
				flag_summ += flag_mass[i];
				if (flag_mass[i] == 1)
					flag_idx = i;
			}
			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] stats_mass;
			delete[] max_L_mass;

			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};

		// классификация по значениям статистики хи-квадрат (с объединением интервалов чобы было больше 5 элементов)
		auto chi_square_stats = [&](double* data, int data_size) {
			int i, j, k, mix_params_amount, buf_intervals_amount;
			double chi_stat, teor_nu, quant_chi;
			bool flag = true;
			int intervals_amount = 30;
			double max_value = find_k_stat(data, data_size, data_size - 1) + 1;
			double min_value = find_k_stat(data, data_size, 0) - 1;
			double len_interval = (max_value - min_value) / intervals_amount;
			if (mixture_type == NORMAL)
				mix_params_amount = 2;
			if (mixture_type == RAYLEIGH)
				mix_params_amount = 1;
			if (mixture_type == LOGNORMAL)
				mix_params_amount = 2;
			double* nu_i = new double[intervals_amount];
			double* interval_bounds = new double[intervals_amount];
			double* nu_i_bounds = new double[intervals_amount];
			double *max_L_mass = new double[hyp_cl_amount];
			int* flag_mass = new int[hyp_cl_amount];

			for (i = 0; i < intervals_amount; ++i)
				nu_i[i] = 0;
			for (i = 0; i < data_size; ++i) {
				k = 0;
				flag = true;
				while (flag) {
					if (((len_interval * k + min_value) <= data[i]) && ((len_interval * (k + 1) + min_value) > data[i]))
						flag = false;
					else
						k++;
				}
				nu_i[k] = nu_i[k] + 1;
			}
			flag = true;
			int r = data_size;
			int buf_nu = 0;
			buf_intervals_amount = 0;
			for (int j = 0; j < intervals_amount; ++j) {
				//cout << "nu_i[k] " << nu_i[j] << endl;
				if ((nu_i[j] > 5) && flag) {
					if (j < intervals_amount - 1) {
						if ((r - nu_i[j] > 5)) {
							interval_bounds[buf_intervals_amount] = len_interval * (j + 1) + min_value;
							nu_i_bounds[buf_intervals_amount] = nu_i[j];
							r -= nu_i[j];
							buf_intervals_amount++;
						}
						else {
							flag = false;
							buf_nu = nu_i[j];
							r -= nu_i[j];
						}
					}
					else {
						interval_bounds[buf_intervals_amount] = len_interval * (j + 1) + min_value;
						nu_i_bounds[buf_intervals_amount] = nu_i[j];
						r -= nu_i[j];
						buf_intervals_amount++;
					}
				}
				else {
					if (flag) {
						buf_nu = nu_i[j];
						r -= nu_i[j];
						flag = false;
					}
					else {
						buf_nu += nu_i[j];
						r -= nu_i[j];
					}
					if (buf_nu > 5) {
						if (j < intervals_amount - 1) {
							if (r > 5)
							{
								interval_bounds[buf_intervals_amount] = len_interval * (j + 1) + min_value;
								nu_i_bounds[buf_intervals_amount] = buf_nu;
								buf_intervals_amount++;
								buf_nu = 0;
								flag = true;
							}
						}
						else {
							//cout << "j " << j << " buf_nu " << buf_nu << endl;
							interval_bounds[buf_intervals_amount] = len_interval * (j + 1) + min_value;
							nu_i_bounds[buf_intervals_amount] = buf_nu;
							buf_intervals_amount++;
							buf_nu = 0;
							flag = true;
						}
					}
				}

			}

			int flag_summ = 0;
			int flag_idx = 0;

			if (buf_intervals_amount > mix_params_amount + 1) {
				//cout << "+" << endl;
				for (i = 0; i < hyp_cl_amount; ++i) {
					chi_stat = 0;

					for (j = 0; j < buf_intervals_amount; ++j) {
						if (j == 0) {
							if (mixture_type == NORMAL)
								teor_nu = cdf(normal(mix_shift[i], mix_scale[i]), interval_bounds[j]);
							if (mixture_type == LOGNORMAL)
								teor_nu = cdf(lognormal(mix_shift[i], mix_scale[i]), interval_bounds[j]);
							if (mixture_type == RAYLEIGH)
								teor_nu = cdf(rayleigh(mix_scale[i]), interval_bounds[j]);
						}
						else {
							if (j != buf_intervals_amount - 1) {
								if (mixture_type == NORMAL)
									teor_nu = cdf(normal(mix_shift[i], mix_scale[i]), interval_bounds[j])
									- cdf(normal(mix_shift[i], mix_scale[i]), interval_bounds[j - 1]);
								if (mixture_type == LOGNORMAL)
									teor_nu = cdf(lognormal(mix_shift[i], mix_scale[i]), interval_bounds[j])
									- cdf(lognormal(mix_shift[i], mix_scale[i]), interval_bounds[j - 1]);
								if (mixture_type == RAYLEIGH)
									teor_nu = cdf(rayleigh(mix_scale[i]), interval_bounds[j])
									- cdf(rayleigh(mix_scale[i]), interval_bounds[j - 1]);
							}
							else {
								if (mixture_type == NORMAL)
									teor_nu = 1 - cdf(normal(mix_shift[i], mix_scale[i]), interval_bounds[j - 1]);
								if (mixture_type == LOGNORMAL)
									teor_nu = 1 - cdf(lognormal(mix_shift[i], mix_scale[i]), interval_bounds[j - 1]);
								if (mixture_type == RAYLEIGH)
									teor_nu = 1 - cdf(rayleigh(mix_scale[i]), interval_bounds[j - 1]);
							}
						}
						//cout << "teor_nu " << teor_nu << endl;
						teor_nu = teor_nu * (data_size);
						chi_stat += (nu_i_bounds[j] - teor_nu)* (nu_i_bounds[j] - teor_nu) / teor_nu;

					}
					quant_chi = quantile(chi_squared(buf_intervals_amount - 1 - mix_params_amount), 0.99);
					//cout << "classs " << i << ": chi_stat - " << chi_stat << " quant_chi = " << quant_chi << endl;
					if (chi_stat < quant_chi)
						flag_mass[i] = 1;
					else
						flag_mass[i] = 0;
				}

				for (i = 0; i < hyp_cl_amount; ++i) {
					flag_summ += flag_mass[i];
					if (flag_mass[i] == 1)
						flag_idx = i;
				}
			}
			else
				flag_summ = 0;

			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] nu_i;
			delete[] interval_bounds;
			delete[] nu_i_bounds;
			delete[] max_L_mass;

			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};

		// классификация по значениям статистики Крамера-мизеса
		auto kramer_mizes_smirnoff_stats = [&](double* data, int data_size,
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			unsigned i, j;
			//int intervals_amount = 60;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr;
			int* flag_mass = new int[hyp_cl_amount];
			double* stats_mass = new double[hyp_cl_amount];
			double* Ui_mass = new double[data_size];
			int flag_summ = 0;
			int flag_idx = 0;
			int unic_amount = 1;
			double last_elem;
			double *max_L_mass = new double[hyp_cl_amount];
			quickSort(data, 0, data_size - 1, data_size / 2);
			
			double dn_bound = 0.4614;
			
			for (i = 0; i < hyp_cl_amount; ++i) {
				F_n_curr = 0;
				buf_d_n = 0;
				
				for (j = 0; j < data_size; ++j) {
					if (mixture_type == LOGNORMAL)
						buf_d_n += pow(cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
					else {
						if (mixture_type == NORMAL)
							buf_d_n += pow(cdf(normal(mix_shift[i], mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						else {
							if (mixture_type == RAYLEIGH)
								buf_d_n += pow(cdf(rayleigh(mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						}
					}
				}
				buf_d_n += 1.0 / double(12 * data_size);
				buf_d_n = (buf_d_n  - 0.4 / double(data_size) + 0.6 / double(data_size*data_size))*(1 + 1.0 / double(data_size));
				//buf_d_n = (buf_d_n - 0.03 / double(data_size) )*(1 + 0.5 / double(data_size));
				//cout << "buf_d_n - " << buf_d_n << endl;
				if (buf_d_n < dn_bound)
					flag_mass[i] = 1;
				else
					flag_mass[i] = 0;
			}
			for (i = 0; i < hyp_cl_amount; ++i) {
				flag_summ += flag_mass[i];
				if (flag_mass[i] == 1)
					flag_idx = i;
			}
			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] stats_mass;
			delete[] max_L_mass;
			delete[] Ui_mass;
			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};
		
		// классификация по значениям статистики Ватсона
		auto watson_stats = [&](double* data, int data_size, 
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			unsigned i, j;
			//int intervals_amount = 60;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr, last_elem;
			int* flag_mass = new int[hyp_cl_amount];
			double* stats_mass = new double[hyp_cl_amount];
			double* Ui_mass = new double[data_size];
			int flag_summ = 0;
			int flag_idx = 0;
			int unic_amount = 1;
			double *max_L_mass = new double[hyp_cl_amount];
			quickSort(data, 0, data_size - 1, data_size / 2);
			
			double dn_bound = 0.186880;

			for (i = 0; i < hyp_cl_amount; ++i) {
				F_n_curr = 0;
				buf_d_n = 0;
				
				for (j = 0; j < data_size; ++j) {
					if (mixture_type == LOGNORMAL) {
						buf_d_n += pow(cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]) - ( (j + 1) - 0.5) / double( data_size), 2);
						F_n_curr += cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]);
					}
					else {
						if (mixture_type == NORMAL)
							buf_d_n += pow(cdf(normal(mix_shift[i], mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						else {
							if (mixture_type == RAYLEIGH)
								buf_d_n += pow(cdf(rayleigh(mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						}
					}
				}
				buf_d_n += 1.0 / double(12 * data_size)- data_size*pow(F_n_curr/double(data_size)-0.5,2);
				//buf_d_n = (buf_d_n / 4.0 - 0.4 / data_size + 0.6 / double(data_size*data_size))*(1 + 1.0 / double(data_size));
				//cout << "buf_d_n - " << buf_d_n << endl;
				if (buf_d_n < dn_bound)
					flag_mass[i] = 1;
				else
					flag_mass[i] = 0;
			}
			for (i = 0; i < hyp_cl_amount; ++i) {
				flag_summ += flag_mass[i];
				if (flag_mass[i] == 1)
					flag_idx = i;
			}
			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] stats_mass;
			delete[] max_L_mass;
			delete[] Ui_mass;
			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};

		// классификация по значениям статистики Андерсона
		auto anderson_stats = [&](double* data, int data_size,
				shared_ptr<double[]> mix_shift, shared_ptr<double[]> mix_scale, int  hyp_cl_amount) {
			int i, j;
			//int intervals_amount = 60;
			int k, buf_intervals_amount;
			bool flag = true;
			double max_d_n, buf_d_n, F_n_curr;
			int* flag_mass = new int[hyp_cl_amount];
			double* stats_mass = new double[hyp_cl_amount];
			double* Ui_mass = new double[data_size];
			int flag_summ = 0;
			int flag_idx = 0;
			int unic_amount = 1;
			double last_elem;
			double *max_L_mass = new double[hyp_cl_amount];
			quickSort(data, 0, data_size - 1, data_size / 2);

			double dn_bound = 2.4924;

			for (i = 0; i < hyp_cl_amount; ++i) {
				F_n_curr = 0;
				buf_d_n = 0;

				for (j = 0; j < data_size; ++j) {
					if (mixture_type == LOGNORMAL) {
						buf_d_n += log(cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]))*((2*(j+1)-1)/(double(2*data_size)))
							+ (-(2*(j + 1)-1) / double(2*data_size)+1)*log(1.0-cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]));
						//F_n_curr += cdf(lognormal(mix_shift[i], mix_scale[i]), data[j]);
					}
					else {
						if (mixture_type == NORMAL)
							buf_d_n += pow(cdf(normal(mix_shift[i], mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						else {
							if (mixture_type == RAYLEIGH)
								buf_d_n += pow(cdf(rayleigh(mix_scale[i]), data[j]) - (2 * (j + 1) - 1) / double(2 * data_size), 2);
						}
					}
				}
				buf_d_n = -buf_d_n*2- data_size;
				//buf_d_n = (buf_d_n / 4.0 - 0.4 / data_size + 0.6 / double(data_size*data_size))*(1 + 1.0 / double(data_size));
				//cout << "buf_d_n - " << buf_d_n << endl;
				if (buf_d_n < dn_bound)
					flag_mass[i] = 1;
				else
					flag_mass[i] = 0;
			}
			for (i = 0; i < hyp_cl_amount; ++i) {
				flag_summ += flag_mass[i];
				if (flag_mass[i] == 1)
					flag_idx = i;
			}
			if (flag_summ > 1) {
				for (int i = 0; i < hyp_cl_amount; ++i)
					L_max_calculation(data, data_size, i, mix_shift[i], mix_scale[i], max_L_mass);

				int beg_idx = 0;
				int max_idx = 0;

				for (int m = beg_idx; m < hyp_cl_amount; ++m) {
					if (max_L_mass[max_idx] < max_L_mass[m])
						max_idx = m;
				}
				flag_summ = 1;
				flag_idx = max_idx;
			}
			delete[] flag_mass;
			delete[] stats_mass;
			delete[] max_L_mass;
			delete[] Ui_mass;
			if (flag_summ == 1)
				return flag_idx;
			else return -1;
		};
		
		// копирование одномерного массива в двумерный
		auto copy_in_one_mass = [&](double* image_one_mass, int x_c, int y_c, int x_l, int y_l) {
			int idx = 0;
			int i, j;
			for (i = x_c; i < x_c + x_l; ++i) {
				for (j = y_c; j < y_c + y_l; ++j) {
					image_one_mass[idx] = raw_image[img_l_x*i +j];
					idx++;
				}
			}
		};
		x_l = length_;
		#pragma omp for
		for (int i = 0; i < iters_i_x + 1; ++i) {
			if (i == iters_i_x)
				x_l = length_ + (img_l_x - (iters_i_x * half_step + length_));
			y_l = length_;
			for ( int j = 0; j < iters_i_y + 1; ++j) {
				
				if (j == iters_i_y)
					y_l = length_ + (img_l_y - (iters_i_y * half_step + length_));
				
				copy_in_one_mass(buf_img, i*(half_step), j*(half_step), x_l, y_l);

				//idx_class = chi_square_stats(buf_img,  x_l*y_l);

				
				idx_class = kolmogorov_stats2(buf_img, x_l*y_l, mix_shift, mix_scale, loc_hyp_cl_amount);
				//idx_class = kramer_mizes_smirnoff_stats(buf_img, x_l*y_l, mix_shift, mix_scale, loc_hyp_cl_amount);
				//idx_class = watson_stats(buf_img, x_l*y_l, mix_shift, mix_scale, loc_hyp_cl_amount);
				//idx_class = anderson_stats(buf_img, x_l*y_l, mix_shift, mix_scale, loc_hyp_cl_amount);
				//idx_class = max_likehood_stats(buf_img, x_l*y_l, mix_shift, mix_scale, hyp_cl_amount);
				//cout << "idx_class " << idx_class <<"  i*(length_/2), j*(length_/2), "<< i * (half_step)<< " "<< j*(half_step)<< endl;
				if (idx_class != -1) {
					#pragma omp critical
					{
						for (l = i * (half_step); l < i*(half_step)+x_l; ++l)
							for (m = j * (half_step); m < j*(half_step)+y_l; ++m)
								class_flag[l][m] = idx_class + 1;
					}
				}
			}
		}
		delete[] buf_img;
	}
	auto end1 = std::chrono::steady_clock::now();
	auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
	cout << "elapsed_ms1  " << elapsed_ms1.count() << endl;
}

//BIC

void mixture_handler::BIC() {
	double summ = 0;
	double big_summ = 0;
	long double f_i_j, pix_buf;
	unsigned i, j, idx_i, idx_j;
	

	for (i = 0; i < img_l_x*img_l_y; ++i) {
		idx_i = i / img_l_y;
		idx_j = i % img_l_y;
		summ = 0;
		pix_buf = raw_image[idx_i*img_l_x+idx_j];
		for (j = 0; j < hyp_cl_amount; ++j)
			summ += mix_weight[j] * (1 / (mix_scale[j] * sqrt(2 * pi)))*exp(-((pix_buf
				- mix_shift[j])*(pix_buf - mix_shift[j])) / (2.0 * mix_scale[j] * mix_scale[j]));

		big_summ += log(summ);
	}
	unsigned count_n_z = 0;
	/*for (j = 0; j < hyp_cl_amount; ++j)
		if (mix_weight[j] != 0)
			++count_n_z;*/
	count_n_z = hyp_cl_amount;
	bic_value = -2 * big_summ + log(img_l_x*img_l_y)*(3 * count_n_z - 1);

	cout << "BIC:  " << bic_value << "     " << big_summ << "\n" << "\n" << endl;
	
}

// загрузка результат классификации в файл 

void mixture_handler::printInformation_to_image() {
	CImage result;
	result.Create(img_l_y, img_l_x, 24);
	int cl_idx_color_g = 0;
	int cl_idx_color_r = 0;
	int cl_idx_color_b = 0;
	// задаем цвет пикселя
	for (int i = 0; i < img_l_x; i++) {
		for (int j = 0; j < img_l_y; j++) {
			cl_idx_color_g = (255/hyp_cl_amount)* class_flag[img_l_x - i - 1][j];
			cl_idx_color_b = (255 / (hyp_cl_amount/2))*(class_flag[img_l_x - i - 1][j] % 3);
			cl_idx_color_r = (255 / (hyp_cl_amount / 3))*(class_flag[img_l_x - i - 1][j] % 3);
			result.SetPixelRGB(j, i, cl_idx_color_r, cl_idx_color_g, cl_idx_color_b);
		}
	}

	CString _name = L"D:\\classification_image_.jpg";
	/*LPCTSTR file_name = LPCTSTR(_name.c_str());
	cout << file_name << endl;*/
	result.Save(_name);

}

void mixture_handler::detect_result_by_mask(string filename) {
	std::ofstream vmdelet_out;     //создаем поток 
	vmdelet_out.open(filename, std::ios::app);
	CImage mask_image;
	double amount_cl_pix = 0;
	unsigned curr_class;
	double amount_true_pix = 0;
	int y_len, x_len, i, j, k;
	cout << "percentage of correctly classified pixels:" << endl;
	for (k = 0; k < hyp_cl_amount; ++k) {
		if (img_mask_list[k] != "\"\"")
		{
			mask_image.Load(img_mask_list[k].c_str());
			amount_cl_pix = 0;
			amount_true_pix = 0;
			for (i = 0; i < img_l_x; i++) {
				for (j = 0; j < img_l_y; j++) {
					curr_class = (unsigned(GetBValue(mask_image.GetPixel(j, i ))) / 255)*(k + 1);
					if (curr_class != 0) {
						amount_cl_pix++;
						//cout << i << " " << j << class_flag[img_l_x - i - 1][j] << curr_class + 1 << endl;
						if (class_flag[img_l_x - i - 1][j] ==( curr_class+1))
							//if (class_flag[i][j] ==( curr_class + 1))
								//if (class_flag[j][i] == (curr_class + 1))
							amount_true_pix++;
					}

				}
			}
			mask_image.Detach();
			cout << "class " << k + 1 << ", "<< curr_class<< ": " << amount_true_pix / amount_cl_pix << endl;
			
			// открываем файл для записи 
			
			vmdelet_out << "class " << k + 1 << ", " << curr_class << ": " << amount_true_pix / amount_cl_pix << "\n"; // сама запись
			
		}	
	}
	vmdelet_out.close();   // закрываем файл
}

// вывод информации о полученном результате классификации

void mixture_handler::printInformation() {
	unsigned i, j;
	cout << "finded model:" << "\n";
	/*cout << " mixture components: " << hyp_cl_amount << "\n";
	cout << "EM mix_shift values:" << "\n";
	for (i = 0; i < hyp_cl_amount; i++)
		cout << mix_shift[i] << "  ";
	cout << endl;
	cout << endl;
	cout << "EM mix_scale values:" << "\n";
	for (i = 0; i < hyp_cl_amount; i++)
		cout << mix_scale[i] << "  ";
	cout << endl;*/

	out.open(split_mix_filename);
	for (i = 0; i < img_l_x; i++) {
		for (j = 0; j < img_l_y; j++)
			out << class_flag[i][j] << " ";
		out << std::endl;
	}
	out.close();
}

//поиск медианы

double mixture_handler::find_med(double* window, int wind_size) {
	int med_index = (wind_size) / 2 - 1;
	bool flag = true;
	int left = 0;
	int right = wind_size - 1;
	if (med_index >= 0) {
		while (flag) {
			//cout << "ddd";
			std::pair<int, int> result = partition(window, left, right, med_index);
			if (result.first< med_index && result.second > med_index) {
				flag = false;

			}
			else {
				if (result.first > med_index)
					right = result.first;
				else {
					if (result.second < med_index)
						left = result.second;
				}
			}
		}

		return 	window[med_index];
	}
	else
		return -1;
}

//поиск к-той порядковой статистики

double mixture_handler::find_k_stat(double * data, int wind_size, int k_stat) {
	bool flag = true;
	int  left = 0;
	int  right = wind_size - 1;

	while (flag) {
		std::pair<int, int> result = partition(data, left, right, k_stat);
		if (result.first< k_stat && result.second > k_stat)
			flag = false;
		else {
			if (result.first > k_stat)
				right = result.first;
			else {
				if (result.second < k_stat)
					left = result.second;
			}
		}
	}

	return 	data[k_stat];
}

// быстрая сортировка

void  mixture_handler::quickSort(double * data,  int l, int r, int pivot_index) {
	double v , temp;
	int i ,	j ,	p ,	q ;
	while (l < r) {
		v = data[r];
		i = l;
		j = r - 1;
		p = l - 1;
		q = r;
		while (i <= j) {
			while ((data[i] < v)&&(i<r))
				i++;
			while ((data[j] > v)&&(j>0))
				j--;
			if (i >= j)
				break;

			temp = data[i];
			data[i] = data[j];
			data[j] =temp;
			//swap(data[i], data[j]);
			if (data[i] == v) {
				p++;
				temp = data[p];
				data[p] = data[i];
				data[i] = temp;
				//swap(data[i], data[p]);
			}
			i++;
			if (data[j] == v) {
				q--;
				temp = data[q];
				data[q] = data[j];
				data[j] = temp;
				//swap(data[q], data[j]);
			}
			j--;
		}
		temp = data[i];
		data[i] = data[r];
		data[r] = temp;
		//swap(data[i], data[r]);
		j = i - 1;
		i++;
		for (int k = l; k <= p; k++, j--) {
			temp = data[k];
			data[k] = data[j];
			data[j] = temp;
			//swap(data[k], data[j]);
		}

		for (int k = r - 1; k >= q; k--, i++) {
			temp = data[k];
			data[k] = data[i];
			data[i] = temp;
			//swap(data[i], data[k]);
		}

		if ((j - l )<( r - i)){
			quickSort(data, l, j, 0);
			l = i;
		}
		else {
			quickSort(data, i, r, 0);
			r = j;
		}
	}
}


std::pair<int, int> mixture_handler::partition(double* mass, int left, int right, int  ind_pivot)
{
	double v = mass[ind_pivot];
	//	double buf = mass[left];
	//double v = mass[right];

	double temp = mass[right];
	mass[right] = v;
	mass[ind_pivot] = temp;
	int i = left;
	int j = right - 1;
	int p = left - 1;
	int q = right;
	while (i <= j) {
		while (mass[i] < v)
			i++;
		while (mass[j] > v)
			j--;
		if (i >= j)
			break;
		
		temp = mass[i];
		mass[i] = mass[j];
		mass[j] = temp;
		if (mass[i] == v) {
			p++;
			temp = mass[p];
			mass[p] = mass[i];
			mass[i] = temp;
		}
		i++;
		if (mass[j] == v) {
			q--;
			temp = mass[q];
			mass[q] = mass[j];
			mass[j] = temp;
		}
		j--;
	}
	temp = mass[i];
	mass[i] = mass[right];
	mass[right] = temp;
	j = i - 1;
	i++;
	for (int k = left; k <= p; k++, j--) {
		temp = mass[k];
		mass[k] = mass[j];
		mass[j] = temp;
	}
		
	for (int k = right - 1; k >= q; k--, i++) {
		temp = mass[k];
		mass[k] = mass[i];
		mass[i] = temp;
	}
		
	return  std::pair<int, int>(j , i);

}

// вызов скрипта на python для отрисовки результатов

void mixture_handler::draw_graphics() {
	string cmd = "echo python  C:\\Users\\anastasya\\PycharmProjects\\untitled5\\mixture_vizualization.py " + gen_mix_filename +
		" " + split_mix_filename +
		" | %windir%\\system32\\cmd.exe \"/K\" C:\\Users\\anastasya\\Anaconda3\\Scripts\\activate.bat  ";
	system(cmd.c_str());
}

// деструктор


