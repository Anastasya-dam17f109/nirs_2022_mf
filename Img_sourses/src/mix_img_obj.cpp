
#include <fstream> // ���������� �����
#include "mix_img_obj.h"

boost::random::mt19937 generator_{ static_cast<std::uint32_t>(time(0)) };

// �����������, ����������� �������� �������������� ���������� �����������

mix_img_obj::mix_img_obj(int img_size, mix_type mix_t, int amount_targets, int classes) {
	image_len_x = img_size * amount_targets;
	image_len_y = img_size * amount_targets;
	alloc_layer_mmr();
	mixture_type = mix_t;
	amount_trg = amount_targets;
	class_amount = classes;
    img_generator();
	img_accumulation();
	print_results();
}

// �����������, ����������� �������� �������������� ���������� ����������� ��� �������� ��� �� �����

mix_img_obj::mix_img_obj(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes) {
	image_len_x = img_size * amount_targets;
	image_len_y = img_size * amount_targets;
	alloc_layer_mmr();
	mixture_type = mix_t;
	amount_trg = amount_targets;
	class_amount = classes;
	img_generator_from_file(file_name);
	img_accumulation();
	print_results();
}


// �����������, ����������� ������������ �������� ��� ,����������� �� �����

mix_img_obj::mix_img_obj(string file_name) {
	ifstream load_params;
	int files_amount = 0;
	int n_buf;
	int mask_amount = 0;
	int stat_mask_amount = 0;
	string buf;
	int input_mix_num = 0;
	
	load_params.open(file_name);
	load_params >> files_amount;
	load_params >> class_amount;
	load_params >> mask_amount;
	load_params >> stat_mask_amount;
	//load_params >> input_mix_num;
	//getline(load_params, buf);
	//mixture_type = static_cast<mix_type>(input_mix_num);;

	class_list         = unique_ptr<wstring[]>(new wstring[class_amount]);
	load_mask_list     = make_unique<wstring[]>(mask_amount);

	mask_list    = shared_ptr<wstring[]>(new wstring[stat_mask_amount]);
	re_mix_shift = shared_ptr<double[]>(new double[class_amount]);
	re_mix_scale = shared_ptr<double[]>(new double[class_amount]);

	load_params >> buf;
	if (buf == "lognormal")
		mixture_type = LOGNORMAL;
	load_params >> buf;
	filename_load_image = wstring(buf.begin(), buf.end());
	
	for (unsigned i = 0; i < class_amount; ++i) {
		load_params >> buf;
		class_list[i] = wstring(buf.begin(), buf.end());
	}

	for (int i = 0; i < mask_amount; ++i) {
		load_params >> buf;
		load_mask_list[i] = wstring(buf.begin(), buf.end());
	}

	for (int i = 0; i < mask_amount; ++i) {
		load_params >> n_buf;
		load_mask_list_idx.push_back(n_buf);
	}

	for (int i = 0; i < stat_mask_amount; ++i) {
		load_params >> buf;
		mask_list[i] = wstring(buf.begin(), buf.end());
	}
	load_params.close();
	cout << "im_here" << endl;
	load_from_bitmap();
	
	img_accumulation();
	print_results();
}

// �������� ����������� �� �����

void mix_img_obj::load_from_bitmap() {
	int y_len, x_len,  j, k, buf_amount=0;
	unsigned i;
	bool mask_flag = false;
	CImage image;
	CImage image_mask;

	image.Load(filename_load_image.c_str());
	y_len = image.GetHeight();
	x_len = image.GetWidth();
	// ������ ����������� ������ 2^n*2^n
	if (y_len > x_len)
		image_len_x = x_len;
	else
		image_len_x = y_len;
	i = 1;
	while (i < image_len_x)
		i *= 2;
	i /= 2;
	image_len_x = i;
	image_len_y = i;
	alloc_layer_mmr();
	cout << "im_here" << endl;
	for (i = 0; i < image_len_x; i++) 
		for (j = 0; j < image_len_y; j++) 
			layer_mx_img[layer_amount-1][image_len_x - i - 1][j] = double(GetGValue(image.GetPixel(j, i)));

	image.Detach();
	for (k = 0; k < class_amount; ++k) {
		image.Load(class_list[k].c_str());
		re_mix_shift[k] = 0;
		re_mix_scale[k] = 0;
		buf_amount = 0;
		y_len = image.GetHeight();
		x_len = image.GetWidth();
		auto pred = find(load_mask_list_idx.begin(), load_mask_list_idx.end(), k);
		if (pred != load_mask_list_idx.end()) {
			mask_flag = true;
			cout << distance(load_mask_list_idx.begin(), pred) << endl;
			image_mask.Load(load_mask_list[distance(load_mask_list_idx.begin(), pred)].c_str());
		}
		else
			mask_flag = false;
		
		for (i = 0; i < y_len; i++)
			for (j = 0; j < x_len; j++) {
				if (mask_flag) {
					re_mix_shift[k] += double(GetGValue(image.GetPixel(j, i)))*(1 - int(GetGValue(image_mask.GetPixel(j, i))) / 255);
					buf_amount += (1 - int(GetGValue(image_mask.GetPixel(j, i))) / 255);
				}
				else {
					re_mix_shift[k] += double(GetGValue(image.GetPixel(j, i)));
					buf_amount = y_len * x_len;
				}
			}
		
		re_mix_shift[k] /= double(buf_amount);
		//cout << "buf_amount " << buf_amount << " " << y_len * x_len << endl;
		for (i = 0; i < y_len; i++)
			for (j = 0; j < x_len; j++) {
				if (mask_flag) 
					re_mix_scale[k] += pow((double((GetGValue(image.GetPixel(j, i)))) - re_mix_shift[k])*(1 - int(GetGValue(image_mask.GetPixel(j, i))) / 255), 2);
				else
					re_mix_scale[k] += pow(double(GetGValue(image.GetPixel(j, i))) - re_mix_shift[k], 2);
			}
		
		re_mix_scale[k] = re_mix_scale[k] / double(buf_amount);
		re_mix_scale[k] = sqrt(log(re_mix_scale[k] / (re_mix_shift[k] * re_mix_shift[k]) + 1.0));
		re_mix_shift[k] = log(re_mix_shift[k] / exp(re_mix_scale[k] * re_mix_scale[k] / 2.0));
		image.Detach();
	}
}

//�������� �������� � ��������� �����������

void mix_img_obj::img_generator() {
	boost::random::normal_distribution <> dist_norm_bcg{ 128, 37.0 };
	boost::random::normal_distribution <> dist_norm_trg{ 240, 1.5 };
	boost::random::uniform_01 <> dist_rel;

	unsigned i, j, k, l, bright_step, amount_brigh_trg, t_coord_x, t_coord_y;
	unsigned mix_number = 1;
	double sred;
	unsigned * targ_size;
	double * targ_bright;

	auto dist_gen_bcg = [&]() {
		if (mixture_type == NORMAL)
			return dist_norm_bcg(generator_);

		else {
			if (mixture_type == RAYLEIGH)
				return sqrt(-2 * pow(20, 2.0) *log(1 - dist_rel(generator_)));
		}
	};

	auto dist_gen_trg = [&](unsigned i) {
		if (mixture_type == NORMAL) {
			if (class_amount+1 > 1)
				dist_norm_trg.param(boost::random::normal_distribution <>::param_type(targ_bright[i], re_mix_scale[i+1]));
			return dist_norm_trg(generator_);
		}
		else {
			if (mixture_type == RAYLEIGH)
				return sqrt(-2 * pow(40, 2.0) *log(1 - dist_rel(generator_)));
		}
	};

	re_mix_shift = shared_ptr<double[]>(new double[class_amount]);
	re_mix_scale = shared_ptr<double[]>(new double[class_amount]);
	for (i = 0; i < class_amount; ++i) {
		re_mix_shift[i] = 0.0;
		re_mix_scale[i] = 0.0;
	}

	targs = shared_ptr<target[]>(new target[amount_trg*amount_trg]);
	targ_size = new unsigned[amount_trg];
	targ_bright = new double[amount_trg];

	
	for (i = 0; i < image_len_x; i++) {
		for (j = 0; j < image_len_x; j++)
			layer_mx_img[layer_amount - 1][i][j] = dist_gen_bcg();
	}

	sred = mean(layer_mx_img[layer_amount - 1]);
	bright_step = (255 - sred - 40) / class_amount;
	amount_brigh_trg = amount_trg / class_amount;
	if (amount_brigh_trg == 0)
		amount_brigh_trg = 1;

	for (i = 0; i < amount_trg; i++) {
		targ_size[i] = min_targ_size + i * 2;
		targ_bright[i] = sred + 40 + (unsigned(i / amount_brigh_trg) + 1) * bright_step;
	}
	re_mix_shift[0] = 128.0;
	re_mix_scale[0] = 37.0;
	re_mix_shift[1] = targ_bright[0];
	re_mix_scale[1] = 1.5;
	if (mixture_type == NORMAL) {
		re_mix_shift[0] = 128.0;
		re_mix_scale[0] = 87.0;
		//re_mix_shift[1] = re_targ_shift;
		re_mix_scale[1] = 80.0;
	}
	else {
		if (mixture_type == RAYLEIGH) {
			re_mix_shift[0] = 0;
			re_mix_scale[0] = 20;
			re_mix_shift[1] = 0;
			//re_targ_shift = 40.0;
			re_mix_scale[1] = 40.0;
		}
	}
	targ_bright[0] = 180;
	for (i = 0; i < amount_trg; i++) {
		if (mixture_type == NORMAL) {
			if (i > 0 && targ_bright[i] != targ_bright[i - 1]) {
				mix_number++;
				re_mix_shift[mix_number] = targ_bright[i];
				//re_mix_scale[mix_number] = 1.5*mix_number;
			}
		}
		for (j = 0; j < amount_trg; j++) {
			t_coord_x = i * backg_size + backg_size / 2 - 1 - targ_size[j] / 2;

			t_coord_y = j * backg_size + backg_size / 2 - 1 - targ_size[j] / 2;

			if (mixture_type == NORMAL)
				targs[i*amount_trg + j].brightness = targ_bright[i];
			else
				targs[i*amount_trg + j].brightness = 30;
			targs[i*amount_trg + j].size = targ_size[j];
			targs[i*amount_trg + j].x = t_coord_x;
			targs[i*amount_trg + j].y = t_coord_y;

			for (k = 0; k < targ_size[j]; k++) {
				for (l = 0; l < targ_size[j]; l++) 
					layer_mx_img[layer_amount - 1][t_coord_x + k][t_coord_y + l] = dist_gen_trg(i);
			}
		}
	}

	for (i = 0; i < amount_trg*amount_trg; i++) {
		for (j = 1; j < class_amount; j++)
			if (targs[i].brightness == re_mix_shift[j])
				targs[i].mix_type = j + 1;
	}
	if (mixture_type == NORMAL) {
		re_mix_shift[0] = 128.0;
		re_mix_scale[0] = 87.0;
		re_mix_shift[1] = 180;
		re_mix_scale[1] = 80.0;
	}
	if (mixture_type == RAYLEIGH) {
		re_mix_shift[0] = 20.0;
		re_mix_scale[0] = 20.0;
		re_mix_shift[1] = 40;
		re_mix_scale[1] = 40;
	}
	delete[] targ_size;
	delete[] targ_bright;
}

//�������� ��������� �������� � ��������� �����������

void  mix_img_obj::img_generator_from_file(string file_name) {

	unsigned i, j, k, l, bright_step, amount_brigh_trg, t_coord_x, t_coord_y;
	unsigned mix_number = 1;
	double sred;
	unsigned * targ_size;
	double * targ_bright;

	re_mix_shift = shared_ptr<double[]>(new double[class_amount]);
	re_mix_scale = shared_ptr<double[]>(new double[class_amount]);
	for (i = 0; i < class_amount; ++i) {
		re_mix_shift[i] = 0.0;
		re_mix_scale[i] = 0.0;
	}

	targs = shared_ptr<target[]>(new target[amount_trg*amount_trg]);
	targ_size = new unsigned[amount_trg];
	targ_bright = new double[amount_trg];

	ifstream in;
	in.open(file_name);
	for (i = 0; i < image_len_x; i++) {
		for (j = 0; j < image_len_x; j++)
			in>>layer_mx_img[layer_amount - 1][i][j] ;
	}

	in.close();
	
	
	amount_brigh_trg = amount_trg / class_amount;
	if (amount_brigh_trg == 0)
		amount_brigh_trg = 1;

	for (i = 0; i < amount_trg; i++) {
		targ_size[i] = min_targ_size + i * 2;
		
	}
	re_mix_shift[0] = 128.0;
	re_mix_scale[0] = 37.0;
	re_mix_shift[1] = targ_bright[0];
	re_mix_scale[1] = 1.5;
	if (mixture_type == NORMAL) {
		re_mix_shift[0] = 128.0;
		targ_bright[0] = 180.0;
		re_mix_scale[0] = 87.0;
		re_mix_shift[1] = 180;
		re_mix_scale[1] = 80.0;
	}
	else {
		if (mixture_type == RAYLEIGH) {
			re_mix_shift[0] = 0;
			re_mix_scale[0] = 20;
			re_mix_shift[1] = 0;
			re_mix_scale[1] = 40.0;
		}
	}
	targ_bright[0] = 180;
	for (i = 0; i < amount_trg; i++) {
		if (mixture_type == NORMAL) {
			if (i > 0 && targ_bright[i] != targ_bright[i - 1]) {
				mix_number++;
				re_mix_shift[mix_number] = targ_bright[i];
				re_mix_scale[mix_number] = 1.5*mix_number;
			}
		}
		for (j = 0; j < amount_trg; j++) {
			t_coord_x = i * backg_size + backg_size / 2 - 1 - targ_size[j] / 2;

			t_coord_y = j * backg_size + backg_size / 2 - 1 - targ_size[j] / 2;

			if (mixture_type == NORMAL)
				targs[i*amount_trg + j].brightness = targ_bright[i];
			else
				targs[i*amount_trg + j].brightness = 30;
			targs[i*amount_trg + j].size = targ_size[j];
			targs[i*amount_trg + j].x = t_coord_x;
			targs[i*amount_trg + j].y = t_coord_y;
		}
	}

	for (i = 0; i < amount_trg*amount_trg; i++) {
		for (j = 1; j < class_amount; j++)
			if (targs[i].brightness == re_mix_shift[j])
				targs[i].mix_type = j + 1;
	}
	if (mixture_type == RAYLEIGH) {
		re_mix_shift[0] = 20.0;
		re_mix_scale[0] = 20.0;
		re_mix_shift[1] = 40;
		re_mix_scale[1] = 40;
	}
	if (mixture_type == NORMAL) {
		re_mix_shift[0] = 128.0;
		targ_bright[0] = 180.0;
		re_mix_scale[0] = 87.0;
		re_mix_shift[1] = 180;
		re_mix_scale[1] = 80.0;
	}
	delete[] targ_size;
	delete[] targ_bright;
}

// ��������� ������ ��� �������������� �����������

void  mix_img_obj::alloc_layer_mmr(){
	int i = 2;
	while (i < image_len_x) {
		i *= 2;
		layer_amount++;
	}

	layer_size = shared_ptr <int[]>(new int[layer_amount]);
	layer_mx_img = new double**[layer_amount];
		//shared_ptr<shared_ptr<shared_ptr<double[]>[]>[]>(new shared_ptr<shared_ptr<double[]>[]>[layer_amount]);

	for (int i = 0; i < layer_amount; ++i) {
		layer_size[i] = pow(2, i + 1);
		layer_mx_img[i] = new double*[layer_size[i]];
			//shared_ptr<shared_ptr<double[]>[]>(new shared_ptr<double[]>[layer_size[i]]);
		for (int j = 0; j < layer_size[i]; j++)
			layer_mx_img[i][j] = new double[layer_size[i]];
			//shared_ptr<double[]>(new double[layer_size[i]]);
	}
}

// ���������� ������������� ����� ������������
// ����������� �������� ���������� �������� ������� ����������

void  mix_img_obj::img_accumulation() {
	for (int i = layer_amount - 2; i > -1; i--) 
		for (int k = 0; k < layer_size[i]; ++k) 
			for (int l = 0; l < layer_size[i]; ++l) {
				layer_mx_img[i][k][l] = (layer_mx_img[i + 1][2 * k][2 * l] + layer_mx_img[i + 1][2 * k + 1][2 * l]
					+ layer_mx_img[ i + 1][2 * k][2 * l + 1] + layer_mx_img[i + 1][2 * k + 1][2 * l + 1]) / 4.0;
				//cout << k << " " << l << endl;
			}
}

//����� �������� �� �����������

void  mix_img_obj::print_results() {
	unsigned i, j;
	out.open(filename_gen_image);
	for (i = 0; i < image_len_x; i++) {
		for (j = 0; j < image_len_y; j++)
			out << layer_mx_img[layer_amount -1][i][j] << " ";
		out << std::endl;
	}
	out.close();
	cout << " generated image params:" << "\n";
	//cout << "mixture type: " << mixture_type << "\n";
	cout << "size: " << image_len_x << " " << image_len_y << "\n";
	cout << " mix components amount: " << class_amount  << "\n";
	cout << "re_mix_shift values:" << endl;
	for (i = 0; i < class_amount; i++)
		cout << re_mix_shift[i] << "  ";
	cout << "\n";

	cout << "re_mix_scale values:" << "\n";
	for (i = 0; i < class_amount; i++)
		cout << re_mix_scale[i] << "  ";
	cout << endl;
}

//���������� �������� ���������������

int mix_img_obj::mean(double** data) {
	double result = 0;
	for (int k = 0; k < image_len_x; k++) {
		for (int l = 0; l < image_len_y; l++)
			result += data[k][l];
	}
	return int(result / (image_len_x*image_len_y));
}
