#include "pch.h"
#include "quad_tree_handler.h"


quad_tree_handler::quad_tree_handler(shared_ptr<mix_img_obj> image){
	m_image = image;
	int i = 2;
	while (i < m_image->get_image_len().first) {
		i *= 2;
		layer_amount++;
	}

	layer_size  = shared_ptr <int[]>(new int[layer_amount]);
	layer_order = shared_ptr<Hilbert_curve[]>(new Hilbert_curve[layer_amount]);
	layer       = shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]>(new shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[layer_amount]);

	for (int i = 0; i < layer_amount; ++i) {
		layer_size[i] = pow(2, i + 1);
		layer[i] = shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>(new shared_ptr<shared_ptr<Node>[]>[layer_size[i]]);
		for (int j = 0; j < layer_size[i]; j++)
			layer[i][j] = shared_ptr<shared_ptr<Node>[]>(new shared_ptr<Node>[layer_size[i]]);
		layer_order[i] = Hilbert_curve(layer_size[i]);
	}
	m_root = shared_ptr<Node>(new Node);
	m_root->parent = nullptr;
	m_root->l_corner = pair<int, int>(0, 0);
	m_root->width = 1;
	build_quad_tree(m_root, 0);
	for (int i = 0; i < layer_amount; ++i) {
		
		for (int j = 0; j < layer_size[i]; j++) {
			for (int k = 0; k < layer_size[i]; k++)
				cout << "(" << layer[i][j][k]->l_corner.first << " " << layer[i][j][k]->l_corner.second << ") ";
			cout << endl;
		}
	}
}


void quad_tree_handler::build_quad_tree(shared_ptr<Node> elem, int n_layer) {
	if (n_layer == layer_amount-1) {
		for (int i = 0; i < 4; ++i) {
			elem->m_children[i] = shared_ptr<Node>(new Node);
			elem->m_children[i]->parent = elem;
			elem->m_children[i]->width = elem->m_children[i]->parent->width;
		}
		
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				layer[n_layer][elem->l_corner.first*2 + i][elem->l_corner.second*2 + j] = elem->m_children[j + 2 * i];
				elem->m_children[j + 2 * i]->l_corner
					= pair<int, int>(elem->l_corner.first * 2 + i, elem->l_corner.second * 2 + j);
			}
		}
				
	}
	else {
		for (int i = 0; i < 4; ++i) {
			elem->m_children[i] = shared_ptr<Node>(new Node);
			elem->m_children[i]->parent = elem;
			elem->m_children[i]->width = elem->m_children[i]->parent->width * 2;
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

quad_tree_handler::~quad_tree_handler()
{
}
