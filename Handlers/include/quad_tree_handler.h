#pragma once
#include "initial_prob_img.h"
#include "Zig_zag_curve.h"
#include "Hilbert_curve.h"
#include "helix_curve.h"


struct Node
{
	shared_ptr<Node> parent;
	array<shared_ptr<Node>, 4> m_children;
	pair<int, int> l_corner;
	

	// ��� ��������� ��� ���� ������� 
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
	double**** init_prob_img;

    unsigned ** class_flag;
    pair<int,int> l_coner;
	shared_ptr<Node>  m_root;
	unsigned layer_amount = 1;
	const unsigned layer_ord_amount = 10;
	shared_ptr<shared_ptr<shared_ptr<shared_ptr<Node>[]>[]>[]> layer;
	shared_ptr <shared_ptr<shared_ptr<Basic_curve>[]>[]> layer_order;
	shared_ptr <int[]>          layer_size;
	shared_ptr<shared_ptr<double[]>[]> p_xs_xs1;
	shared_ptr<shared_ptr<double[]>[]> p_xs_layer;
	double theta = 0.7;
	unsigned class_amount = 1;

	shared_ptr <int[]>          pix_cl_amount;
	shared_ptr <unsigned[]>     idx_max;
	shared_ptr <double[]>       buf_max;
	
	string filename_gen_image = "D:\\generated_image.txt";
	string filename_split_image = "D:\\splitted_image.txt";
public:
	
    quad_tree_handler(shared_ptr<initial_prob_img> image, unsigned** cnt);
    quad_tree_handler(shared_ptr<initial_prob_img> image, int size, unsigned** cnt);
	void build_quad_tree(shared_ptr<Node> elem, int n_layer);
    void set_dest_cnt(unsigned ** cl_fl_ptr, int size);
	void set_probabilities(int i_idx, int j_idx);
	void p_xs_matrix_generator();
	void p_xs_layer_generator();
	void bottom_up_pass();
	void up_down_pass();
	void split_image_by_summ();
	void split_image_by_vote();
	void create_splitted_img();
	void draw_graphics();
    void clear_mem();
    ~quad_tree_handler(){}

};

