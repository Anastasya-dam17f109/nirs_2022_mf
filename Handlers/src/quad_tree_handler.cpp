//#include "pch.h"
#include "quad_tree_handler.h"


quad_tree_handler::quad_tree_handler(shared_ptr<basic_prob_img> image, unsigned** cnt, int _h_tree){
	m_image = image;
	init_prob_img   =  m_image->get_image();
	init_layer_idx  = m_image->get_init_layer_idx();
	init_layer_size = m_image->get_init_layer_size();
	class_amount = m_image->get_class_amount();
    set_dest_cnt(cnt, m_image->get_image_len().first);
	h_tree =  _h_tree;
	int i = 2;
	while (i < m_image->get_image_len().first) {
		i *= 2;
		layer_amount++;
	}

	layer_size  = shared_ptr <int[]>(new int[layer_amount]);
	layer_idx = shared_ptr <int[]>(new int[layer_amount]);
	layer_order = shared_ptr < shared_ptr<shared_ptr<Basic_curve>[]>[]>(new  shared_ptr<shared_ptr<Basic_curve>[]>[layer_ord_amount]);
	for(int i = 0; i< layer_ord_amount; ++i)
		layer_order[i] = shared_ptr<shared_ptr<Basic_curve>[]>(new shared_ptr<Basic_curve>[layer_amount]);
	
	//layer       = shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]>(new shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[layer_amount]);

	p_xs_xs1   = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[class_amount]);
	p_xs_layer = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[layer_amount]);

	int summ = 0;
	for (int i = 0; i < layer_amount; ++i) {
		layer_size[i] = pow(2, i + 1);
		layer_idx[i]  = summ;
		summ += layer_size[i] * layer_size[i];
		//layer[i] = shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>(new shared_ptr<shared_ptr<Node>[]>[layer_size[i]]);
		//for (int j = 0; j < layer_size[i]; j++)
			//layer[i][j] = shared_ptr<shared_ptr<Node>[]>(new shared_ptr<Node>[layer_size[i]]);

		// получение порядков - последовательности пикселей - на каждом слое
		
		/*for (int j = 0; j < 4; ++j) {
			layer_order[j][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], j, false));
			layer_order[j][i]->get_points_for_curve();
		}*/
		layer_order[0][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 0, false));
		layer_order[0][i]->get_points_for_curve();
		layer_order[1][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 0, true));
		layer_order[1][i]->get_points_for_curve();
		
		layer_order[2][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 2, false));
		layer_order[2][i]->get_points_for_curve();
		layer_order[3][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 2, true));
		layer_order[3][i]->get_points_for_curve();
	
		// спирали
		layer_order[4][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 2, false));
		layer_order[4][i]->get_points_for_curve();
		layer_order[5][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 0, false));
		layer_order[5][i]->get_points_for_curve();

		layer_order[8][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 1, false));
		layer_order[8][i]->get_points_for_curve();

		layer_order[6][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 0, false));
		layer_order[6][i]->get_points_for_curve();
		layer_order[7][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 0, true));
		layer_order[7][i]->get_points_for_curve();
	
		layer_order[9][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 2, true));
		layer_order[9][i]->get_points_for_curve();
		/*layer_order[10][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 2, false));
		layer_order[10][i]->get_points_for_curve();*/

		
		/*layer_order[6][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],3, false));
		layer_order[6][i]->get_points_for_curve();
		
		layer_order[7][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],3,true));
		layer_order[7][i]->get_points_for_curve();
		*/

		
		
		/*layer_order[10][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 1, false));
		layer_order[10][i]->get_points_for_curve();
		layer_order[11][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 3, false));
		layer_order[11][i]->get_points_for_curve();*/
		

		

		/*layer_order[11][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],1,true));
		layer_order[11][i]->get_points_for_curve();
		*/

		/*for (int j = 6; j < 6+4; ++j) {
			layer_order[j][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], j-6));
			layer_order[j][i]->get_points_for_curve();
			layer_order[j][i]->reverse_curve();
		}*/
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
	layer = new Node[summ];
	m_root = shared_ptr<Node>(new Node);
	m_root->parent = nullptr;
	m_root->l_corner = pair<int, int>(0, 0);

	build_quad_tree();
	p_xs_matrix_generator();
	p_xs_layer_generator();
}

//

quad_tree_handler::quad_tree_handler(shared_ptr<basic_prob_img> image, int size, unsigned** cnt, int _h_tree) {
	m_image = image;
	h_tree = _h_tree;
	init_prob_img = m_image->get_image();
	init_layer_idx = m_image->get_init_layer_idx();
	init_layer_size = m_image->get_init_layer_size();
	class_amount = m_image->get_class_amount();
	pix_cl_amount = shared_ptr <int[]>(new int[class_amount]);
	idx_max = shared_ptr <unsigned[]>(new unsigned[layer_ord_amount]);
	buf_max = shared_ptr <double[]>(new double[layer_ord_amount]);
    set_dest_cnt(cnt,size);
	int i = 2;
	while (i < size) {
		i *= 2;
		layer_amount++;
	}

	layer_size = shared_ptr <int[]>(new int[layer_amount]);
	layer_idx = shared_ptr <int[]>(new int[layer_amount]);
	layer_order = shared_ptr < shared_ptr<shared_ptr<Basic_curve>[]>[]>(new  shared_ptr<shared_ptr<Basic_curve>[]>[layer_ord_amount]);
	for (int i = 0; i < layer_ord_amount; ++i)
		layer_order[i] = shared_ptr<shared_ptr<Basic_curve>[]>(new shared_ptr<Basic_curve>[layer_amount]);

	//layer = shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]>(new shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[layer_amount]);

	p_xs_xs1 = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[class_amount]);
	p_xs_layer = shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[layer_amount]);
	int summ = 0;
	for (int i = 0; i < layer_amount; ++i) {
		layer_size[i] = pow(2, i + 1);
		layer_idx[i] = summ;
		summ += layer_size[i] * layer_size[i];
		//layer[i] = shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>(new shared_ptr<shared_ptr<Node>[]>[layer_size[i]]);
		//for (int j = 0; j < layer_size[i]; j++)
			//layer[i][j] = shared_ptr<shared_ptr<Node>[]>(new shared_ptr<Node>[layer_size[i]]);

		// получение порядков - последовательности пикселей - на каждом слое

		/*for (int j = 0; j < 4; ++j) {
			layer_order[j][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], j, false));
			layer_order[j][i]->get_points_for_curve();
		}*/
		layer_order[0][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 0, false));
		layer_order[0][i]->get_points_for_curve();
		layer_order[1][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 0, true));
		layer_order[1][i]->get_points_for_curve();

		layer_order[2][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 2, false));
		layer_order[2][i]->get_points_for_curve();
		layer_order[3][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], 2, true));
		layer_order[3][i]->get_points_for_curve();

		// спирали
		layer_order[4][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 2, false));
		layer_order[4][i]->get_points_for_curve();
		layer_order[5][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 0, false));
		layer_order[5][i]->get_points_for_curve();

		layer_order[8][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 1, false));
		layer_order[8][i]->get_points_for_curve();

		layer_order[6][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 0, false));
		layer_order[6][i]->get_points_for_curve();
		layer_order[7][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 0, true));
		layer_order[7][i]->get_points_for_curve();

		layer_order[9][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 2, true));
		layer_order[9][i]->get_points_for_curve();
		/*layer_order[10][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i], 2, false));
		layer_order[10][i]->get_points_for_curve();*/


		/*layer_order[6][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],3, false));
		layer_order[6][i]->get_points_for_curve();

		layer_order[7][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],3,true));
		layer_order[7][i]->get_points_for_curve();
		*/

		/*layer_order[10][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 1, false));
		layer_order[10][i]->get_points_for_curve();
		layer_order[11][i] = shared_ptr<Helix_curve>(new Helix_curve(layer_size[i], 3, false));
		layer_order[11][i]->get_points_for_curve();*/

		/*layer_order[11][i] = shared_ptr<Hilbert_curve>(new Hilbert_curve(layer_size[i],1,true));
		layer_order[11][i]->get_points_for_curve();
		*/

		/*for (int j = 6; j < 6+4; ++j) {
			layer_order[j][i] = shared_ptr<Zig_zag_curve>(new Zig_zag_curve(layer_size[i], j-6));
			layer_order[j][i]->get_points_for_curve();
			layer_order[j][i]->reverse_curve();
		}*/
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

	layer = new Node[summ];
	m_root = shared_ptr<Node>(new Node);
	m_root->parent = nullptr;
	m_root->l_corner = pair<int, int>(0, 0);

	build_quad_tree();
	p_xs_matrix_generator();
	p_xs_layer_generator();
}

//назначение контейнера, в который будет писаться результат классификации

void quad_tree_handler::set_dest_cnt(unsigned ** cl_fl_ptr, int size){
    if(cl_fl_ptr == nullptr){
        class_flag = new unsigned*[size];
        for (int i = 0; i < m_image->get_image_len().first; ++i)
            class_flag[i] = new unsigned[size];
    }
    else
        class_flag =  cl_fl_ptr;

}

//установка значений начальных вероятностей
//передается координата блока!

void quad_tree_handler::set_probabilities(int i_idx, int j_idx) {
    l_coner.first = i_idx;
    l_coner.second = j_idx;
	
}

// построение структуры квадродерева, причем так, чтобы пиксели группировались в слои

//void quad_tree_handler::build_quad_tree(shared_ptr<Node> elem, int n_layer) {
//	if (n_layer == layer_amount-1) {
//		for (int i = 0; i < 4; ++i) {
//			elem->m_children[i] = shared_ptr<Node>(new Node);
//			elem->m_children[i]->parent = elem;
//			
//			elem->m_children[i]->p_xs_ds = unique_ptr<long double[]>(new long double[class_amount]);
//			elem->m_children[i]->p_xs_cs_ds = unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]>(new unique_ptr<unique_ptr<long double[]>[]>[class_amount]);
//			elem->m_children[i]->p_xs_Y = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[layer_ord_amount]);
//            for (unsigned j = 0; j < layer_ord_amount; ++j)
//				elem->m_children[i]->p_xs_Y[j] = unique_ptr<long double[]>(new long double[class_amount]);
//            for (unsigned j = 0; j < class_amount; ++j) {
//				elem->m_children[i]->p_xs_ds[j] = 0;
//                for(unsigned k = 0; k < layer_ord_amount; ++k)
//					elem->m_children[i]->p_xs_Y[k][j] = 0;
//				elem->m_children[i]->p_xs_cs_ds[j] = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[class_amount]);
//                for (unsigned k = 0; k < class_amount; ++k) {
//					elem->m_children[i]->p_xs_cs_ds[j][k] = 0;
//					elem->m_children[i]->p_xs_cs_ds[j][k] = unique_ptr<long double[]>(new long double[class_amount]);
//				}
//			}
//		}
//		
//		for (int i = 0; i < 2; ++i) {
//			for (int j = 0; j < 2; ++j) {
//				layer[layer_idx[n_layer]+
//					layer_size[n_layer]*(elem->l_corner.first * 2 + i)+
//					elem->l_corner.second * 2 + j
//				]
//				//layer[n_layer][elem->l_corner.first*2 + i][elem->l_corner.second*2 + j]
//					= elem->m_children[j + 2 * i];
//				elem->m_children[j + 2 * i]->l_corner
//					= pair<int, int>(elem->l_corner.first * 2 + i, elem->l_corner.second * 2 + j);
//				//for (int k = 0; k < class_amount; ++k) {
//					/*elem->m_children[j + 2 * i]->p_xs_ys[k]
//						= m_image->get_image()[n_layer][elem->m_children[j + 2 * i]->l_corner.first][elem->m_children[j + 2 * i]->l_corner.second][k];*/
//					//elem->m_children[j + 2 * i]->p_xs_ds[k]
//						//= elem->m_children[j + 2 * i]->p_xs_ys[k];
//				//}
//			}
//		}
//				
//	}
//	else {
//		for (int i = 0; i < 4; ++i) {
//			elem->m_children[i] = shared_ptr<Node>(new Node);
//			elem->m_children[i]->parent = elem;
//			
//			elem->m_children[i]->p_xs_ds = unique_ptr<long double[]>(new long double[class_amount]);
//			elem->m_children[i]->p_xs_cs_ds = unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]>(new unique_ptr<unique_ptr<long double[]>[]>[class_amount]);
//			elem->m_children[i]->p_xs_Y = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[layer_ord_amount]);
//			for (int j = 0; j < layer_ord_amount; ++j)
//				elem->m_children[i]->p_xs_Y[j] = unique_ptr<long double[]>(new long double[class_amount]);
//			for (int j = 0; j < class_amount; ++j) {
//				
//				elem->m_children[i]->p_xs_ds[j] = 0;
//				for (int k = 0; k < layer_ord_amount; ++k)
//					elem->m_children[i]->p_xs_Y[k][j] = 0;
//				elem->m_children[i]->p_xs_cs_ds[j] = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[class_amount]);
//				for (int k = 0; k < class_amount; ++k) {
//					elem->m_children[i]->p_xs_cs_ds[j][k] = 0;
//					elem->m_children[i]->p_xs_cs_ds[j][k] = unique_ptr<long double[]>(new long double[class_amount]);
//				}
//			}
//			//elem->m_children[i]->width = elem->m_children[i]->parent->width * 2;
//		}
//
//		for (int i = 0; i < 2; ++i) {
//			for (int j = 0; j < 2; ++j) {
//				layer[layer_idx[n_layer] +
//					layer_size[n_layer] * (elem->l_corner.first * 2 + i) +
//					elem->l_corner.second * 2 + j
//				]
//				
//					= elem->m_children[j + 2 * i];
//				elem->m_children[j + 2 * i]->l_corner
//					= pair<int, int>(elem->l_corner.first* 2 + i, elem->l_corner.second* 2 + j);
//				/*for (int k = 0; k < class_amount; ++k) {
//					elem->m_children[j + 2 * i]->p_xs_ys[k]
//						= m_image->get_image()[n_layer][elem->m_children[j + 2 * i]->l_corner.first][elem->m_children[j + 2 * i]->l_corner.second][k];
//					
//				}*/
//			}
//		}
//		for (int i = 0; i < 4; ++i)
//			build_quad_tree(elem->m_children[i], n_layer + 1);
//	}
//}

void quad_tree_handler::build_quad_tree() {
	for (int i = 0; i < layer_amount; ++i) {
		for (int j = 0; j < layer_size[i]; ++j) {
			for (int k = 0; k < layer_size[i]; ++k) {
				layer[layer_idx[i] +layer_size[i]*j +k].p_xs_ds = unique_ptr<long double[]>(new long double[class_amount]);
				layer[layer_idx[i] + layer_size[i] * j + k].observed = unique_ptr<bool[]>(new bool[layer_ord_amount]);
				layer[layer_idx[i] + layer_size[i] * j + k].p_xs_cs_ds = unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]>
						(new unique_ptr<unique_ptr<long double[]>[]>[class_amount]);
				layer[layer_idx[i] + layer_size[i] * j + k].p_xs_Y = unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[layer_ord_amount]);
				for (int l = 0; l < layer_ord_amount; ++l) {
					layer[layer_idx[i] + layer_size[i] * j + k].p_xs_Y[l] = unique_ptr<long double[]>(new long double[class_amount]);
					layer[layer_idx[i] + layer_size[i] * j + k].observed[l] = false;
				}
				for (int l = 0; l < class_amount; ++l) {
					layer[layer_idx[i] + layer_size[i] * j + k].p_xs_cs_ds[l] =
						unique_ptr<unique_ptr<long double[]>[]>(new unique_ptr<long double[]>[class_amount]);
					for (int t = 0; t < class_amount; ++t)
						layer[layer_idx[i] + layer_size[i] * j + k].p_xs_cs_ds[l][t] = unique_ptr<long double[]>(new long double[class_amount]);
				}


				if (i != layer_amount-1) {
					layer[layer_idx[i] + layer_size[i] * j + k].m_children[0] = &layer[layer_idx[i+1] + layer_size[i+1] * 2*j + 2*k];
					layer[layer_idx[i] + layer_size[i] * j + k].m_children[1] = &layer[layer_idx[i + 1] + layer_size[i + 1] * 2 * j + 2 * k+1];
					layer[layer_idx[i] + layer_size[i] * j + k].m_children[2] = &layer[layer_idx[i + 1] + layer_size[i + 1] * (2 * j+1)+2* k];
					layer[layer_idx[i] + layer_size[i] * j + k].m_children[3] = &layer[layer_idx[i + 1] + layer_size[i + 1] * (2 * j + 1) +2 * k + 1];
					for (int l = 0; l < 4; ++l)
						layer[layer_idx[i] + layer_size[i] * j + k].m_children[l]->parent = &layer[layer_idx[i] + layer_size[i] * j + k];
				}
			}
		}
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
// реализуем формулы из zerubia
// необходимое условие дя работы - нормировка вероятностей

void quad_tree_handler::bottom_up_pass() {
	long double prod_buf, sum_buf;
	double summ;
	int layer_offset = m_image->get_layer_amount() - layer_amount;
	int i, j, k, l, t, r;
	for ( i = layer_amount - 1; i > -1; i--) {
		// обработка листового слоя : p_xs_ds=p_xs_ys
		//  p_xs_ys уже заданы в initial_prob_img, они тут не вычисляются
		
		if (i != (layer_amount - 1)) {
			// вычисление p_xs_ds на узлах-сучках
			
			for (j = 0; j < layer_size[i]; j++) {
				for (k = 0; k < layer_size[i]; k++) {
					summ = 0;
					for (l = 0; l < class_amount; l++) {
						prod_buf = 1;
						sum_buf = 0;
						/*for (t = 0; t < 4; ++t) {
							sum_buf = 0;
							for ( r = 0; r < class_amount; r++) 
								sum_buf += p_xs_xs1[r][l] * layer[i][j][k]->m_children[t]->p_xs_ds[r] / p_xs_layer[i][r];
							prod_buf *= sum_buf;*/
							for (r = 0; r < class_amount; r++)
								sum_buf += p_xs_xs1[r][l] * 
								layer[layer_idx[i+1] + 2*j * layer_size[i+1] + 2*k].p_xs_ds[r] / p_xs_layer[i][r];
							prod_buf *= sum_buf;
							sum_buf = 0;
							for (r = 0; r < class_amount; r++)
								sum_buf += p_xs_xs1[r][l] * 
								layer[ layer_idx[i+1] + (2*j+1) * layer_size[i+1] +2* k ].p_xs_ds[r] / p_xs_layer[i][r];
							prod_buf *= sum_buf;
							sum_buf = 0;
							for (r = 0; r < class_amount; r++)
								sum_buf += p_xs_xs1[r][l] * layer[layer_idx[i+1] +2* j * layer_size[i+1] + 2*k+1]
								.p_xs_ds[r] / p_xs_layer[i][r];
							prod_buf *= sum_buf;
							sum_buf = 0;
							for (r = 0; r < class_amount; r++)
								sum_buf += p_xs_xs1[r][l] * layer[layer_idx[i + 1] + (2 * j+1) * layer_size[i + 1] + 2 * k + 1]
								.p_xs_ds[r] / p_xs_layer[i][r];
							prod_buf *= sum_buf;
						//}
						layer[layer_idx[i]+j*layer_size[i]+k].p_xs_ds[l] = 
							init_prob_img[init_layer_idx[layer_offset + i]+
							(l_coner.first*layer_size[i] + j)* init_layer_size[layer_offset + i] * class_amount +
							(l_coner.second*layer_size[i] + k)* class_amount +l] * prod_buf;
						/*cout << init_prob_img[init_layer_idx[layer_offset + i] +
							(l_coner.first*layer_size[i] + j)* init_layer_size[layer_offset + i] * class_amount +
							(l_coner.second*layer_size[i] + k)* class_amount + l] << endl;*/
							/*layer[i][j][k]->p_xs_ys[l] * prod_buf;*/
						summ += layer[layer_idx[i] + j * layer_size[i] + k].p_xs_ds[l];
					}
					for (l = 0; l < class_amount; l++)
					{

						layer[layer_idx[i] + j * layer_size[i] + k].p_xs_ds[l] /= summ;
						
					}
					for (l = 0; l < layer_ord_amount; l++)
						layer[layer_idx[i] + j * layer_size[i] + k].observed[l] = false;
				}
			}
		}
		else {
			for (j = 0; j < layer_size[i]; j++)
				for (k = 0; k < layer_size[i]; k++) {
					for (l = 0; l < class_amount; l++)
						//layer[i][j][k]->p_xs_ds[l] = layer[i][j][k]->p_xs_ys[l] ;
						layer[layer_idx[i] + j * layer_size[i] + k].p_xs_ds[l] =
						init_prob_img[init_layer_idx[layer_offset + i] +
						(l_coner.first*layer_size[i] + j)* init_layer_size[layer_offset + i] * class_amount +
						(l_coner.second*layer_size[i] + k)* class_amount + l];
					for (l = 0; l < layer_ord_amount; l++)
						layer[layer_idx[i] + j * layer_size[i] + k].observed[l] = false;
				}
		}
		// p_xs_ds вычислены целиком все
		// вычисление p_xs_ds_cs
		// начинаем с самого начала - с листьев
		for (j = 0; j < layer_size[i]; j++) {
			for (k = 0; k < layer_size[i]; k++) {
				summ = 0;
				for (l = 0; l < class_amount; l++) {
					for (t = 0; t < class_amount; t++) {
						for ( r = 0; r < class_amount; r++) {
							if (i != 0)
								layer[layer_idx[i] + j * layer_size[i] + k].p_xs_cs_ds[l][t][r] 
									= layer[layer_idx[i] + j * layer_size[i] + k].p_xs_ds[l] * p_xs_xs1[l][t]
									* p_xs_xs1[l][r] * p_xs_layer[i][r] * p_xs_layer[i-1][t] / (p_xs_layer[i][l] * p_xs_layer[i][l]);
							else
								layer[layer_idx[i] + j * layer_size[i] + k].p_xs_cs_ds[l][t][r] 
										= layer[layer_idx[i] + j * layer_size[i] + k].p_xs_ds[l]
									* p_xs_xs1[l][r] * p_xs_layer[i][r] / (p_xs_layer[i][l]* p_xs_layer[i][l]);
							summ += layer[layer_idx[i] + j * layer_size[i] + k].p_xs_cs_ds[l][t][r];
						}
					}
				}
                for ( l = 0; l < class_amount; l++)
                    for (t = 0; t < class_amount; t++)
                        for (r = 0; r < class_amount; r++)
							layer[layer_idx[i] + j * layer_size[i] + k].p_xs_cs_ds[l][t][r] /= summ;
				for (l = 0; l < layer_ord_amount; l++)
					layer[layer_idx[i] + j * layer_size[i] + k].observed[l] = false;
			}
		}
	}
}

// проход сверху-вних по квадродереву
// реализуются формулы из zerubia 
// дополнение - все вероятности надо нормировать! иначе переполнение типа
// на самом слое вводится отношение порядка в виде марковской цепи
// тут  реализована только лишь одна гильбертова кривая
// также бльшой вопрос по начальную точку на нулевом слое - как задать инициирующую вероятность

void quad_tree_handler::up_down_pass() {
	Point buf, buf1;
	double sum;
	int l, i, j, k, t, r;
    for (l = 0; l < layer_ord_amount; ++l) {
        for (i = 0; i < layer_amount; i++) {
			// обработка начального слоя - отличие в том, что у этих пикселей нет родителей
			if (i == 0) {
				for (j = 0; j < layer_size[i] * layer_size[i]; j++) {
					buf = layer_order[l][i]->get_points()[j];
					sum = 0;
					
					
					if (j != 0) {
						buf1 = layer_order[l][i]->get_points()[j - 1];
						sum = 0;
						//layer[i][buf1.x][buf1.y]->p_xs_Y[0][0] = 0;
                        for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
                            for (t = 0; t < class_amount; ++t)
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][0][t]
                                        * layer[layer_idx[i] + buf1.x * layer_size[i] + buf1.y].p_xs_Y[l][t];
							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];

						}
						for (k = 0; k < class_amount; ++k)
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
						
					}
					else {

						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
							for (t = 0; t < class_amount; ++t)
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] 
								+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][0][t];

							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
						}
						for (k = 0; k < class_amount; ++k)
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
					}
				}
				/*if (l == 0) {
					cout << "layer: " << i << endl;
					for (int p = 0; p < layer_size[i] * layer_size[i]; ++p) {
						cout << "(" << p / layer_size[i] << "," << p % layer_size[i] << ")" << ": ";
						for (k = 0; k < class_amount; ++k) {

							cout << layer[layer_idx[i] + p].p_xs_Y[l][k] << " ";
						}
						cout << endl;
					}
				}*/
			}
			// обработка промежуточных слоев -  у этих пикселей родители есть, поэтому идет  отличие в формулах
			else {
                for (j = 0; j < layer_size[i] * layer_size[i]; j++) {
					buf = layer_order[l][i]->get_points()[j];
					sum = 0;
					
					
					if (j != 0) {
						buf1 = layer_order[l][i]->get_points()[j - 1];
						sum = 0;
                        //layer[i][buf1.x][buf1.y]->p_xs_Y[0][0] = 0;
                        for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
                            for (t = 0; t < class_amount; ++t)
                                for (r = 0; r < class_amount; ++r)
									layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] 
									+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][r][t]
                                            * layer[layer_idx[i] + buf1.x * layer_size[i] + buf1.y].p_xs_Y[l][t] * layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][r];
							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
						}
                        for ( k = 0; k < class_amount; ++k)
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
					}
					else {

						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
							for (t = 0; t < class_amount; ++t)
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k]
								+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][t][0]
								* layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][t];
							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
						}
						for (k = 0; k < class_amount; ++k)
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
					}
				}
				/*if (l == 0) {
					cout << "layer: " << i << endl;
					for (int p = 0; p < layer_size[i] * layer_size[i]; ++p) {
						cout << "(" << p / layer_size[i] << "," << p % layer_size[i] << ")" << ": ";
						for (k = 0; k < class_amount; ++k) {

							cout << layer[layer_idx[i] + p].p_xs_Y[l][k] << " ";
						}
						cout << endl;
					}
				}*/
			}
		}
	}
}

// проход сверху-вних по квадродереву
// реализуются формулы из zerubia 
// дополнение - все вероятности надо нормировать! иначе переполнение типа
// на самом слое вводится отношение порядка в виде марковской цепи
// тут  реализована только лишь одна гильбертова кривая
// также бльшой вопрос по начальную точку на нулевом слое - как задать инициирующую вероятность
// добавление c однородными площадками

void quad_tree_handler::up_down_pass_V2() {
	Point buf, buf1;
	double sum;
	double max_dev;
	
	int l, i, j, k, t, r;
	for (l = 0; l < layer_ord_amount; ++l) {
		for (i = layer_amount- h_tree; i < layer_amount; i++) {
			// обработка начального слоя - отличие в том, что у этих пикселей нет родителей
			if (i == 0) {
				for (j = 0; j < layer_size[i] * layer_size[i]; j++) {
					buf = layer_order[l][i]->get_points()[j];
					sum = 0;


					if (j != 0) {
						buf1 = layer_order[l][i]->get_points()[j - 1];
						sum = 0;
						//layer[i][buf1.x][buf1.y]->p_xs_Y[0][0] = 0;
						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
							for (t = 0; t < class_amount; ++t)
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][0][t]
								* layer[layer_idx[i] + buf1.x * layer_size[i] + buf1.y].p_xs_Y[l][t];
							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];

						}
						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].observed[l] = true;
						}


					}
					else {

						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
							for (t = 0; t < class_amount; ++t)
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k]
								+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][0][t];

							sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
						}
						for (k = 0; k < class_amount; ++k) {
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
							layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].observed[l] = true;
						}
					}
				}
				/*if (l == 0) {
					cout << "layer: " << i << endl;
					for (int p = 0; p < layer_size[i] * layer_size[i]; ++p) {
						cout << "(" << p / layer_size[i] << "," << p % layer_size[i] << ")" << ": ";
						for (k = 0; k < class_amount; ++k) {

							cout << layer[layer_idx[i] + p].p_xs_Y[l][k] << " ";
						}
						cout << endl;
					}
				}*/
			}
			// обработка промежуточных слоев -  у этих пикселей родители есть, поэтому идет  отличие в формулах
			else {
				for (j = 0; j < layer_size[i] * layer_size[i]; j++) {
					buf = layer_order[l][i]->get_points()[j];
					if (!layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].observed[l]) {
						sum = 0;
						max_dev = 0;

						if (j != 0) {
							buf1 = layer_order[l][i]->get_points()[j - 1];
							sum = 0;
							//layer[i][buf1.x][buf1.y]->p_xs_Y[0][0] = 0;
							for (k = 0; k < class_amount; ++k) {
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
								for (t = 0; t < class_amount; ++t)
									for (r = 0; r < class_amount; ++r)
										layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k]
										+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][r][t]
										* layer[layer_idx[i] + buf1.x * layer_size[i] + buf1.y].p_xs_Y[l][t] * layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][r];
								sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
							}
							for (k = 0; k < class_amount; ++k) {
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].observed[l] = true;
								if (abs(layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] -
									layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][k]) > max_dev)
									max_dev = abs(layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] -
										layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][k]);

							}
						}
						else {

							for (k = 0; k < class_amount; ++k) {
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] = 0;
								for (t = 0; t < class_amount; ++t)
									layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k]
									+= layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_cs_ds[k][t][0]
									* layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][t];
								sum += layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
							}
							for (k = 0; k < class_amount; ++k) {
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] /= sum;
								layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].observed[l] = true;
								if (abs(layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] -
									layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][k]) > max_dev)
									max_dev = abs(layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k] -
										layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].parent->p_xs_Y[l][k]);

							}
						}
						
						if (max_dev < accuracy) {
							int counter = 1;
							int len = 1;
							while (i + counter != layer_amount) {
								len *= 2;
								for (int id_i = 0; id_i < len; id_i++) {
									for (int id_j = 0; id_j < len; id_j++) {
										for (k = 0; k < class_amount; ++k)
											layer[layer_idx[i + counter] + (buf.x*len + id_i) * layer_size[i + counter] + buf.y*len + id_j].p_xs_Y[l][k] =
											layer[layer_idx[i] + buf.x * layer_size[i] + buf.y].p_xs_Y[l][k];
										layer[layer_idx[i + counter] + (buf.x*len + id_i) * layer_size[i + counter] + buf.y*len + id_j].observed[l] = true;
									}

								}
								counter++;
							}
						}
					}
				}

				/*if (l == 0) {
					cout << "layer: " << i << endl;
					for (int p = 0; p < layer_size[i] * layer_size[i]; ++p) {
						cout << "(" << p / layer_size[i] << "," << p % layer_size[i] << ")" << ": ";
						for (k = 0; k < class_amount; ++k) {

							cout << layer[layer_idx[i] + p].p_xs_Y[l][k] << " ";
						}
						cout << endl;
					}
				}*/
			}
		}
	}
}


// непосредственная классификация пикселей по значению максимума соответствующей апостериорной вероятности

void quad_tree_handler::split_image_by_summ() {
	double buf_max;
	double buf_prob;
	int idx_max;
    for (int i = 0; i <  layer_size[layer_amount-1]; i++) {
        for (int j = 0; j <  layer_size[layer_amount-1]; j++) {
			buf_max = 0;
			for (int k = 0; k < layer_ord_amount; ++k)
				
				buf_max += layer[layer_idx[layer_amount - 1] + i * layer_size[layer_amount - 1] + j].p_xs_Y[k][0];
				/*buf_max += layer[layer_amount - 1][i][j]->p_xs_ys[0];*/
			
			idx_max = 0;
			for (int l = 0; l < class_amount; l++) {
				//cout << layer[layer_amount - 1][i][j]->p_xs_Y[l] << endl;
				buf_prob = 0;
				for (int k = 0; k < layer_ord_amount; ++k)
					buf_prob += layer[layer_idx[layer_amount - 1] + i * layer_size[layer_amount - 1] + j].p_xs_Y[k][l];
				/*buf_prob += layer[layer_amount - 1][i][j]->p_xs_ys[l];*/
				if (buf_max < buf_prob) {
					buf_max = buf_prob;
					idx_max = l;
				}
			}
			class_flag[l_coner.first*layer_size[layer_amount - 1]+ i][l_coner.second*layer_size[layer_amount - 1]+j] = idx_max + 1;
		}
	}
}

//

void quad_tree_handler::split_image_by_vote() {
	
	
	int weight_idx_max;
    unsigned k;
	int i, j, l;
 
    for (i = l_coner.first*layer_size[layer_amount - 1]; 
			i < l_coner.first*layer_size[layer_amount - 1] + layer_size[layer_amount-1]; i++) {
        for (j = l_coner.second*layer_size[layer_amount - 1]; 
				j < l_coner.second*layer_size[layer_amount - 1] + layer_size[layer_amount-1]; j++) {
			// зануление 
            /*for (k = 0; k < layer_ord_amount; ++k) {
				buf_max[k] = 0;
				idx_max[k] = 0;
			}*/
            for (k = 0; k < class_amount; ++k)
				pix_cl_amount[k] = 0;

            for (k = 0; k < layer_ord_amount; ++k) {
				
				buf_max[k] = layer[layer_idx[layer_amount - 1]
					+ (i - l_coner.first*layer_size[layer_amount - 1]) * layer_size[layer_amount - 1] 
					+ j - l_coner.second*layer_size[layer_amount - 1]].p_xs_Y[k][0];

				idx_max[k] = 0;
                for (l = 0; l < class_amount; l++) {
					/*cout << layer[layer_idx[layer_amount - 1]
						+ (i - l_coner.first*layer_size[layer_amount - 1]) * layer_size[layer_amount - 1]
						+ j - l_coner.second*layer_size[layer_amount - 1]].p_xs_Y[k][l] << endl;*/
					
					
					if (buf_max[k] < layer[layer_idx[layer_amount - 1]
						+ (i - l_coner.first*layer_size[layer_amount - 1]) * layer_size[layer_amount - 1]
						+ j - l_coner.second*layer_size[layer_amount - 1]].p_xs_Y[k][l]) {
						buf_max[k] = layer[layer_idx[layer_amount - 1]
							+ (i - l_coner.first*layer_size[layer_amount - 1]) * layer_size[layer_amount - 1]
							+ j - l_coner.second*layer_size[layer_amount - 1]].p_xs_Y[k][l];
						idx_max[k] = l;
					}
				}
				pix_cl_amount[idx_max[k]] += 1;
			}
			weight_idx_max = 0;
            for (k = 0; k < class_amount; ++k)
				if (pix_cl_amount[weight_idx_max] < pix_cl_amount[k])
					weight_idx_max = k;
			class_flag[i][j] = weight_idx_max + 1;
		}
	}
}


// печать классифицированного изображения  в файл

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

// отрисовка графики - вызов скрипта на python

void quad_tree_handler::draw_graphics() {
	cout << "end 2" << endl;
	string cmd = "echo python  C:\\Users\\anastasya\\PycharmProjects\\untitled5\\mixture_vizualization.py " + filename_gen_image + " " + filename_split_image +
		" | %windir%\\system32\\cmd.exe \"/K\" C:\\Users\\anastasya\\Anaconda3\\Scripts\\activate.bat  ";
	system(cmd.c_str());
}

// очистка памяти в случае если она была выделена

void quad_tree_handler::clear_mem(){
    for(int i = 0; i< layer_size[layer_amount-1]; ++i)
        delete[] class_flag[i];
    delete[] class_flag;
}
