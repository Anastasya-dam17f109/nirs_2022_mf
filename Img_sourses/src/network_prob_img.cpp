#include "network_prob_img.h"


void newtwork_prob_img::gen_prob_img_from_config(string filename, int mode) {
	cout << "creation" << endl;
	m_image = shared_ptr<mix_img_obj>(new mix_img_obj(filename, false, mode));
	m_image_len_x = m_image->get_image_len().first;
	m_image_len_y = m_image->get_image_len().second;
	m_class_amount = m_image->get_class_amount();
	m_layer_amount = m_image->get_layer_amount();
	m_layer_size = m_image->get_layer_size();
	m_layer_idx = m_image->get_layer_idx();
	m_cl_probs = shared_ptr<double[]>(new double[m_class_amount]);

	for (int i = 0; i < m_class_amount; ++i)
		m_cl_probs[i] = 1.0 / m_class_amount;
	cout << "allocation " << endl;
	alloc_layer_mmr();
	//for (int i = 0; i < m_layer_amount; ++i)
		//cout << m_init_layer_idx[i] << " " << m_layer_size[2 * i] << " " << m_layer_size[2 * i + 1] << endl;
	cout << "allocation ok" << endl;
}

//

void newtwork_prob_img::load_probs_from_file(vector<string> probs_data)
{
	int i, j, k, t;
	double buf;
	std::ifstream load_params;
	//здесь предполагаем, что первый файл - это последний слой, нет пропусков и прочего

	for (k = m_layer_amount - 1; k > int(m_layer_amount) - int(probs_data.size()) - 1; --k)
	{
		load_params.open(probs_data[int(m_layer_amount) - 1 - k]);
		load_params >> buf;
		load_params >> buf;

		for (j = 0; j < m_layer_size[2 * k + 1]; j++) {
			for (i = 0; i < m_layer_size[2 * k]; i++) {
				//cout << i << " " << j << endl;
				for (t = 0; t < m_class_amount; t++) {
					load_params >> buf;
					m_prob_img[m_init_layer_idx[k] +
						i * m_layer_size[2 * k + 1] * m_class_amount +
						j * m_class_amount +
						t] = buf;
				}
			}
		}
		load_params.close();
	}
}
