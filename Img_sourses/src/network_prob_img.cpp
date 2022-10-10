//#include "pch.h"
#include "network_prob_img.h"



//

network_prob_img::network_prob_img(string configFilename)
{
	
	
	m_image = shared_ptr<mix_img_obj>(new mix_img_obj(filename, genFlag));

	image_len_x = m_image->get_image_len().first;
	image_len_y = m_image->get_image_len().second;
	class_amount = m_image->get_class_amount();
}




void network_prob_img::alloc_layer_mmr(){

	init_layer_idx = shared_ptr <int[]>(new int [layer_amount]);
	int summ = 0;
	for (int k = 0; k < layer_amount; ++k) {
		init_layer_idx[k] = summ;
		summ += layer_size[k] * layer_size[k] * class_amount;
	}

	init_prob_img = new double [summ];
	
}

