#pragma once
#include "mix_img_obj.h"
#include "Hilbert_curve.h"

struct Node
{
	shared_ptr<Node> parent;
	array<shared_ptr<Node>, 4> m_children;
	pair<int, int> l_corner;
	int width;

	// для обработки как поля маркова 
	unique_ptr<double[]> p_xs_ys;
	unique_ptr<double[]> p_xs_ds;
	unique_ptr<unique_ptr<double[]>[]> p_xs_cs_ds;
	unique_ptr<double[]> p_xs_Y;
	bool isLeaf() const
	{
		return !static_cast<bool>(this->m_children[0]);
	}
	
};

class quad_tree_handler
{
	shared_ptr<mix_img_obj> m_image;

	shared_ptr<Node>  m_root;
	int layer_amount = 1;
	shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]> layer;
	shared_ptr<Hilbert_curve[]> layer_order;
	shared_ptr <int[]>          layer_size;
	//shared_ptr <int[]>          layer_width;
	
public:
	//quad_tree_handler() {};
	quad_tree_handler(shared_ptr<mix_img_obj> image);
	void build_quad_tree(shared_ptr<Node> elem, int n_layer);
	~quad_tree_handler();
};

