#pragma once
#include "mix_img_obj.h"
#include "basic_prob_img.h"
#include "omp.h"

class newtwork_prob_img: public basic_prob_img
{

public:
    newtwork_prob_img(){}

    // конструктор получает имя файла о изображении, а также флаг, нужно или нет самостоятельно генерить вероятности
    virtual void gen_prob_img_from_config(string filename);//, bool flag);
    virtual void load_probs_from_file(list<string> probs_data);

};


