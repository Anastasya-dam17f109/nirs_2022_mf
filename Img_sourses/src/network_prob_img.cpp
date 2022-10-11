#include "network_prob_img.h"


void newtwork_prob_img::gen_prob_img_from_config(string filename) {

    m_image      = shared_ptr<mix_img_obj>(new mix_img_obj(filename, false));
    m_image_len_x  = m_image->get_image_len().first;
    m_image_len_y  = m_image->get_image_len().second;
    m_class_amount = m_image->get_class_amount();
    m_layer_amount = m_image->get_layer_amount();
    m_layer_size   = m_image->get_layer_size();
    m_layer_idx    = m_image->get_layer_idx();
    m_cl_probs     = shared_ptr<double[]>(new double[m_class_amount]);
    for (int i = 0; i < m_class_amount; ++i)
        m_cl_probs[i] = 1.0 / m_class_amount;
    alloc_layer_mmr();
}

//

void newtwork_prob_img::load_probs_from_file(list<string> probs_data)
{
	//здесь предполагаем, что первый файл - это последний слой, нет пропусков и прочего
	for (int k = m_layer_amount - 1; k > m_layer_amount- probs_data.size(); --k)
	{
		std::ifstream load_params;
		
		load_params.open(probs_data[m_layer_amount - 1 -k]);
		double buf;
		load_params >> buf;
		cout << "buf " << buf << endl;
		load_params >> buf;
		cout << "buf " << buf << endl;
		for (int i = 0; i < m_layer_size[k]; i++) {
			for (int j = 0; j < m_layer_size[k]; j++) {
				
				for (int t = 0; t < m_class_amount; t++) {
					load_params >> buf;
					m_prob_img[m_init_layer_idx[k] +
						i * m_layer_size[k] * m_class_amount +
						j * m_class_amount +
						t] = buf;
				}
			}
		}
		load_params.close();
	}
}
