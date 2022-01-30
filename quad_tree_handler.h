#pragma once
#include "initial_prob_img.h"
#include "Zig_zag_curve.h"
#include "Hilbert_curve.h"

struct Node
{
	shared_ptr<Node> parent;
	array<shared_ptr<Node>, 4> m_children;
	pair<int, int> l_corner;
	

	// для обработки как поля маркова 
	shared_ptr<long double[]> p_xs_ys;
	unique_ptr<long double[]> p_xs_ds;
	unique_ptr<unique_ptr<unique_ptr<long double[]>[]>[]> p_xs_cs_ds;
	unique_ptr<unique_ptr<long double[]>[]> p_xs_Y;
	bool isLeaf() const
	{
		return !static_cast<bool>(this->m_children[0]);
	}
	
};

class quad_tree_handler
{
	shared_ptr<initial_prob_img> m_image;
	shared_ptr<shared_ptr<double[]>[]> class_flag;
	shared_ptr<Node>  m_root;
	int layer_amount = 1;
	int layer_ord_amount = 6;
	shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]> layer;
	shared_ptr <shared_ptr<shared_ptr<Basic_curve>[]>[]> layer_order;
	shared_ptr <int[]>          layer_size;
	shared_ptr<shared_ptr<double[]>[]> p_xs_xs1;
	shared_ptr<shared_ptr<double[]>[]> p_xs_layer;
	double theta = 0.9;
	int class_amount = 1;

	
	string filename_gen_image = "D:\\generated_image.txt";
	string filename_split_image = "D:\\splitted_image.txt";
public:
	
	quad_tree_handler(shared_ptr<initial_prob_img> image);
	void build_quad_tree(shared_ptr<Node> elem, int n_layer);
	void p_xs_matrix_generator();
	void p_xs_layer_generator();
	void bottom_up_pass();
	void up_down_pass();
	void split_image();
	void create_splitted_img();
	void draw_graphics();
	~quad_tree_handler();

};

