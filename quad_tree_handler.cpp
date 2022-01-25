#include "pch.h"
#include "quad_tree_handler.h"


quad_tree_handler::quad_tree_handler(shared_ptr<initial_prob_img> image){
	m_image = image;
	class_amount = m_image->get_class_amount();
	class_flag = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[m_image->get_image_len().first]);
	for (int i = 0; i < m_image->get_image_len().first; ++i)
		class_flag[i] = shared_ptr<double[]>(new double[m_image->get_image_len().first]);
	int i = 2;
	while (i < m_image->get_image_len().first) {
		i *= 2;
		layer_amount++;
	}

	layer_size  = shared_ptr <int[]>(new int[layer_amount]);
	layer_order = shared_ptr<Hilbert_curve[]>(new Hilbert_curve[layer_amount]);
	layer       = shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]>(new shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[layer_amount]);

	p_xs_xs1   = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[class_amount]);
	p_xs_layer = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[layer_amount]);

	for (int i = 0; i < layer_amount; ++i) {
		layer_size[i] = pow(2, i + 1);
		layer[i] = shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>(new shared_ptr<shared_ptr<Node>[]>[layer_size[i]]);
		for (int j = 0; j < layer_size[i]; j++)
			layer[i][j] = shared_ptr<shared_ptr<Node>[]>(new shared_ptr<Node>[layer_size[i]]);
		layer_order[i].n_order = layer_size[i];
		layer_order[i].get_points_for_curve();
		p_xs_layer[i] = shared_ptr<double[]>(new double[class_amount]);
		if (i == 0)
			for (int j = 0; j < class_amount; ++j)
				p_xs_layer[i][j] = m_image->get_cl_probs()[j];
		else
			for (int j = 0; j < class_amount; ++j)
				p_xs_layer[i][j] = 0;
	}
	for (int i = 0; i < class_amount; ++i)
		p_xs_xs1[i] = shared_ptr<double[]>(new double[class_amount]);

	m_root = shared_ptr<Node>(new Node);
	m_root->parent = nullptr;
	m_root->l_corner = pair<int, int>(0, 0);

	build_quad_tree(m_root, 0);
	p_xs_matrix_generator();
	p_xs_layer_generator();
	/*for (int i = 0; i < layer_amount; ++i) {
		
		for (int j = 0; j < layer_size[i]; j++) {
			for (int k = 0; k < layer_size[i]; k++)
				cout << "(" << layer[i][j][k]->l_corner.first << " " << layer[i][j][k]->l_corner.second << ") ";
			cout << endl;
		}
	}*/
	/*for (int j = 0; j < layer_size[layer_amount - 1]; j++) {
		for (int k = 0; k < layer_size[layer_amount - 1]; k++)
			cout << "(" << layer[layer_amount-1][j][k]->p_xs_ys[0]- m_image->get_image()[j][k][0]<< ") ";
		cout << endl;
	}*/
}

// построение структуры квадродерева, причем так, чтобы пиксели группировались в слои

void quad_tree_handler::build_quad_tree(shared_ptr<Node> elem, int n_layer) {
	if (n_layer == layer_amount-1) {
		for (int i = 0; i < 4; ++i) {
			elem->m_children[i] = shared_ptr<Node>(new Node);
			elem->m_children[i]->parent = elem;
			elem->m_children[i]->p_xs_ys = shared_ptr<long double[]>(new long double[class_amount]);
			elem->m_children[i]->p_xs_ds = unique_ptr<long double[]>(new long double[class_amount]);
			elem->m_children[i]->p_xs_cs_ds = unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]>(new unique_ptr<unique_ptr<long double[]>[]>[class_amount]);
			elem->m_children[i]->p_xs_Y = unique_ptr<long double[]>(new long double[class_amount]);
			for (int j = 0; j < class_amount; ++j) {
				//elem->m_children[i]->p_xs_ys[j] = 0;
				elem->m_children[i]->p_xs_ds[j] = 0;
				elem->m_children[i]->p_xs_Y[j] = 0;
				elem->m_children[i]->p_xs_cs_ds[j] = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[class_amount]);
				for (int k = 0; k < class_amount; ++k) {
					elem->m_children[i]->p_xs_cs_ds[j][k] = 0;
					elem->m_children[i]->p_xs_cs_ds[j][k] = unique_ptr<long double[]>(new long double[class_amount]);
				}
			}
		}
		
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				layer[n_layer][elem->l_corner.first*2 + i][elem->l_corner.second*2 + j] = elem->m_children[j + 2 * i];
				elem->m_children[j + 2 * i]->l_corner
					= pair<int, int>(elem->l_corner.first * 2 + i, elem->l_corner.second * 2 + j);
				for (int k = 0; k < class_amount; ++k) {
					elem->m_children[j + 2 * i]->p_xs_ys[k]
						= m_image->get_image()[elem->m_children[j + 2 * i]->l_corner.first][elem->m_children[j + 2 * i]->l_corner.second][k];
					elem->m_children[j + 2 * i]->p_xs_ds[k]
						= elem->m_children[j + 2 * i]->p_xs_ys[k];
				}
			}
		}
				
	}
	else {
		for (int i = 0; i < 4; ++i) {
			elem->m_children[i] = shared_ptr<Node>(new Node);
			elem->m_children[i]->parent = elem;
			elem->m_children[i]->p_xs_ys = shared_ptr<long double[]>(new long double[class_amount]);
			elem->m_children[i]->p_xs_ds = unique_ptr<long double[]>(new long double[class_amount]);
			elem->m_children[i]->p_xs_cs_ds = unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]>(new unique_ptr<unique_ptr<long double[]>[]>[class_amount]);
			elem->m_children[i]->p_xs_Y = unique_ptr<long double[]>(new long double[class_amount]);
			for (int j = 0; j < class_amount; ++j) {
				elem->m_children[i]->p_xs_ys[j] = 0;
				elem->m_children[i]->p_xs_ds[j] = 0;
				elem->m_children[i]->p_xs_Y[j] = 0;
				elem->m_children[i]->p_xs_cs_ds[j] = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[class_amount]);
				for (int k = 0; k < class_amount; ++k) {
					elem->m_children[i]->p_xs_cs_ds[j][k] = 0;
					elem->m_children[i]->p_xs_cs_ds[j][k] = unique_ptr<long double[]>(new long double[class_amount]);
				}
			}
			//elem->m_children[i]->width = elem->m_children[i]->parent->width * 2;
		}

		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				layer[n_layer][elem->l_corner.first*2 + i][elem->l_corner.second*2 + j] = elem->m_children[j + 2 * i];
				elem->m_children[j + 2 * i]->l_corner
					= pair<int, int>(elem->l_corner.first* 2 + i, elem->l_corner.second* 2 + j);
			}
		}
		for (int i = 0; i < 4; ++i)
			build_quad_tree(elem->m_children[i], n_layer + 1);
	}
}

// генерация матрицы переходных вероятностей

void quad_tree_handler::p_xs_matrix_generator() {
	for (int i = 0; i < class_amount; ++i)
		for (int j = 0; j < class_amount; ++j) {
			if(i != j)
				p_xs_xs1[i][j] = (1 - theta) / double(class_amount - 1);
			else
				p_xs_xs1[i][j] = theta;
		}
}

// генерация вероятностей появлений каждого класса на каждом слое квадродерева

void quad_tree_handler::p_xs_layer_generator() {
	for (int i = 1; i < layer_amount; ++i) {
		for (int j = 0; j < class_amount; ++j) {
			for (int k = 0; k < class_amount; ++k)
				p_xs_layer[i][j] += p_xs_xs1[j][k] * p_xs_layer[i - 1][k];
			
		}
	}
}

// проход снизу вверх по квадродереву

void quad_tree_handler::bottom_up_pass() {
	for (int i = layer_amount - 1; i > -1; i--) {
		if (i == (layer_amount - 1)) {
			for (int j = 0; j < layer_size[i]; j++) {
				for (int k = 0; k < layer_size[i]; k++)
					for (int l = 0; l < class_amount; l++) {
						layer[i][j][k]->p_xs_ys[l] *= p_xs_layer[i][l];
						layer[i][j][k]->p_xs_ds[l] = layer[i][j][k]->p_xs_ys[l] ;
						//cout << layer[i][j][k]->p_xs_ys[l] << endl;
					}
			}

		}
		else {
			// вычисление p_xs_ys
			for (int j = 0; j < layer_size[i]; j++) {
				for (int k = 0; k < layer_size[i]; k++) {
					double summ = 0;
					for (int l = 0; l < class_amount; l++) {
						layer[i][j][k]->p_xs_ys[l] = 1;
						for (int t = 0; t < 4; ++t)
							layer[i][j][k]->p_xs_ys[l] *= layer[i][j][k]->m_children[t]->p_xs_ys[l] ;
						summ += layer[i][j][k]->p_xs_ys[l];
					}
					//for (int l = 0; l < class_amount; l++) {
						//layer[i][j][k]->p_xs_ys[l] /= summ;
						//cout << i << "  " << layer[i][j][k]->p_xs_ys[l] << endl;
					//}
				}
			}

			// вычисление p_xs_ds
			for (int j = 0; j < layer_size[i]; j++) {
				for (int k = 0; k < layer_size[i]; k++) {
					double summ = 0;
					for (int l = 0; l < class_amount; l++) {
						long double prod_buf = 1;
						long double sum_buf = 0;
						for (int t = 0; t < 4; ++t) {
							sum_buf = 0;
							for (int r = 0; r < class_amount; r++) {
								if (i == 0) {
									cout << p_xs_xs1[r][l] << endl;
									cout << layer[i][j][k]->m_children[t]->p_xs_ds[r] << endl;
									cout << p_xs_layer[i][r] << endl;
								}
								sum_buf += p_xs_xs1[r][l] * layer[i][j][k]->m_children[t]->p_xs_ds[r] / p_xs_layer[i][r];
							}
							prod_buf *= sum_buf;
						}
						layer[i][j][k]->p_xs_ds[l] = layer[i][j][k]->p_xs_ys[l] * prod_buf;
						summ += layer[i][j][k]->p_xs_ds[l];
						cout << i << "  " << layer[i][j][k]->p_xs_ds[l] << endl;
					}
					/*for (int l = 0; l < class_amount; l++)
						layer[i][j][k]->p_xs_ds[l] /= summ;*/
				}
			}
		}
	
	// вычисление p_xs_ds_cs
			for (int j = 0; j < layer_size[i]; j++) {
				for (int k = 0; k < layer_size[i]; k++) {
					double summ = 0;
					for (int l = 0; l < class_amount; l++) {
						for (int t = 0; t < class_amount; t++) {
							for (int r = 0; r < class_amount; r++) {
								if (i != 0)
									layer[i][j][k]->p_xs_cs_ds[l][t][r] = layer[i][j][k]->p_xs_ds[l] * p_xs_xs1[l][t]
									* p_xs_xs1[l][r] * p_xs_layer[i][r] * p_xs_layer[i-1][t] / (p_xs_layer[i][l] * p_xs_layer[i][l]);
								else
									layer[i][j][k]->p_xs_cs_ds[l][t][r] = layer[i][j][k]->p_xs_ds[l] 
									* p_xs_xs1[l][r] * p_xs_layer[i][r] / (p_xs_layer[i][l]* p_xs_layer[i][l]);
								summ += layer[i][j][k]->p_xs_cs_ds[l][t][r];
							}
						}

					}
					/*for (int l = 0; l < class_amount; l++)
						for (int t = 0; t < class_amount; t++)
							for (int r = 0; r < class_amount; r++)
								layer[i][j][k]->p_xs_cs_ds[l][t][r] /= summ;*/
				}
			}



		
	}
}

// проход сверху-вних по квадродереву

void quad_tree_handler::up_down_pass() {
	Point buf, buf1;
	for (int i = 0; i < layer_amount; i++ ){
		if (i == 0) {
			for (int j = 0; j < layer_size[i] * layer_size[i]; j++) {
				buf = layer_order[i].get_points()[j];
				
				if (j == 0) {

					for (int k = 0; k < class_amount; ++k) {
						for (int t = 0; t < class_amount; ++t)
						layer[i][buf.x][buf.y]->p_xs_Y[k] += layer[i][buf.x][buf.y]->p_xs_cs_ds[k][0][t];
						cout << layer[i][buf.x][buf.y]->p_xs_Y[k] << endl;
					}
				}
				else {
					buf1 = layer_order[i].get_points()[j - 1];
					for (int k = 0; k < class_amount; ++k) {
						for (int t = 0; t < class_amount; ++t)
							layer[i][buf.x][buf.y]->p_xs_Y[k] += layer[i][buf.x][buf.y]->p_xs_cs_ds[k][0][t] * layer[i][buf1.x][buf1.y]->p_xs_Y[t];
						cout << layer[i][buf.x][buf.y]->p_xs_Y[k] << endl;
					}
				}
			}
		}
		else {
			for (int j = 0; j < layer_size[i] * layer_size[i]; j++) {
				buf = layer_order[i].get_points()[j];
				
				if (j == 0) {

					for (int k = 0; k < class_amount; ++k) {
						for (int t = 0; t < class_amount; ++t)
							layer[i][buf.x][buf.y]->p_xs_Y[k] += layer[i][buf.x][buf.y]->p_xs_cs_ds[k][t][0] * layer[i][buf.x][buf.y]->parent->p_xs_Y[t];
						cout << layer[i][buf.x][buf.y]->p_xs_Y[k] << endl;
					}
				}
				else {
					buf1 = layer_order[i].get_points()[j - 1];
					for (int k = 0; k < class_amount; ++k) {
						for (int t = 0; t < class_amount; ++t)
							for (int r = 0; r < class_amount; ++r)
								layer[i][buf.x][buf.y]->p_xs_Y[k] += layer[i][buf.x][buf.y]->p_xs_cs_ds[k][r][t]
								* layer[i][buf1.x][buf1.y]->p_xs_Y[t] * layer[i][buf.x][buf.y]->parent->p_xs_Y[r];
						cout << layer[i][buf.x][buf.y]->p_xs_Y[k] << endl;
					}
				}
			}
		}
	}
}

//

void quad_tree_handler::split_image() {
	for (int i =0; i < m_image->get_image_len().first; i++) {
		for (int j = 0; j < m_image->get_image_len().first; j++) {
			
			double buf_max = layer[layer_amount - 1][i][j]->p_xs_Y[0];
			int idx_max = 0;
			for (int l = 0; l < class_amount; l++) {
				//cout << layer[layer_amount - 1][i][j]->p_xs_Y[l] << endl;
				if (buf_max < layer[layer_amount - 1][i][j]->p_xs_Y[l]) {
					buf_max = layer[layer_amount - 1][i][j]->p_xs_Y[l];
					idx_max = l;
				}
			}
			class_flag[i][j] = idx_max + 1;
		}
	}
}

//

void quad_tree_handler::create_splitted_img() {
	ofstream out;
	out.open(filename_split_image);
	for (int i = 0; i < m_image->get_image_len().first; i++) {
		for (int j = 0; j < m_image->get_image_len().first; j++)
			out << class_flag[i][j] << " ";
		out << std::endl;
	}
	out.close();
}


void quad_tree_handler::draw_graphics() {
	cout << "end 2" << endl;
	string cmd = "echo python  C:\\Users\\anastasya\\PycharmProjects\\untitled5\\mixture_vizualization.py " + filename_gen_image + " " + filename_split_image +
		" | %windir%\\system32\\cmd.exe \"/K\" C:\\Users\\anastasya\\Anaconda3\\Scripts\\activate.bat  ";
	system(cmd.c_str());
}
quad_tree_handler::~quad_tree_handler()
{
}
