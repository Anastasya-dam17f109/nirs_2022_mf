#pragma once
#include "basic_prob_img.h"
#include "omp.h"

class initial_prob_img: public basic_prob_img
{
public:
    initial_prob_img(){}
    virtual void gen_prob_img_with_targs(int img_size, mix_type mix_t, int amount_targets, int classes);
    virtual void gen_prob_img_with_targs_from_file(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes);
    virtual void gen_prob_img_from_obj(shared_ptr<mix_img_obj> image);
    virtual void gen_prob_img_from_config(string filename, int mode);

    void     generate_init_probs_mixtures();
	void     generate_init_probs_mix_opMP_V2();
	void     generate_init_probs_max_apost_opMP_V2();
};

