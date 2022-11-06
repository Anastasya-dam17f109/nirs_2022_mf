#pragma once
#include "mix_img_obj.h"
#include "omp.h"
#include <vector>

class basic_prob_img
{
public:
    shared_ptr<mix_img_obj> m_image;
	// размеры слоев изображений
    shared_ptr <int[]>      m_layer_size;
    shared_ptr <int[]>      m_layer_idx;
	// индексы элементов начала слоев в одномерном массиве вероятностей
    shared_ptr <int[]>      m_init_layer_idx;
    int m_layer_amount = 1;
	// одномерный массив вероятностей
    double* m_prob_img;

    shared_ptr<double[]> mix_shift;
    shared_ptr<double[]> mix_scale;
    shared_ptr<double[]> m_cl_probs;
    unsigned m_image_len_x = 32;
    unsigned m_image_len_y = 32;
    unsigned m_class_amount = 1;
    int m_init_window_size = 8;

    basic_prob_img(){}
	// обработка смесями изображения с автоматически добавляемыми целями (с его генерацией)
    virtual void gen_prob_img_with_targs(int img_size, mix_type mix_t, int amount_targets, int classes){}
	// обработка смесями изображения с автоматически добавляемыми целями (с его загрузкой из файла)
    virtual void gen_prob_img_with_targs_from_file(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes){}
	// обработка смесями изображения , заложенного в объект 
    virtual void gen_prob_img_from_obj(shared_ptr<mix_img_obj> image){}

	// построение изображениия через днные,  заложенные в файле конфигурации
    virtual void gen_prob_img_from_config(string filename){}
    virtual void load_probs_from_file(vector<string> probs_data){}
	// выделение памяти под массив вероятностей
    void  alloc_layer_mmr()
    {
        m_init_layer_idx = shared_ptr <int[]>(new int [m_layer_amount]);
        int summ = 0;
		cout << "m_init_layer_idx[k] "<< endl;
        for (int k = 0; k < m_layer_amount; ++k) {
            m_init_layer_idx[k] = summ;
			cout << m_init_layer_idx[k] << endl;
            summ += m_layer_size[2*k] * m_layer_size[2*k+1] * m_class_amount;
        }
        m_prob_img = new double [summ];
    }

    int                     get_class_amount()    {return m_class_amount;}
    std::pair<int, int>     get_image_len()       {return std::pair<int, int>(m_image_len_x, m_image_len_y); }
    double*                 get_image()           {return m_prob_img;}
    shared_ptr<mix_img_obj> get_m_image()         {return m_image;}
    shared_ptr<double[]>    get_cl_probs()        {return m_cl_probs;}
    int                     get_layer_amount()    {return m_layer_amount;}
    shared_ptr <int[]>      get_init_layer_idx () {return m_init_layer_idx;}
    shared_ptr <int[]>      get_init_layer_size() {return m_layer_size;}
    virtual ~basic_prob_img()
    {
        delete [] m_prob_img;
    }
};

