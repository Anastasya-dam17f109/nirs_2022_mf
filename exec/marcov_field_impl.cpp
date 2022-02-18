// marcov_field_impl.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include "quad_tree_handler.h"
#include <boost/filesystem.hpp>
#include "mixture_handler.h"
using namespace boost::filesystem;
int main() {
	//for (int order = 1; order < 3; order++) {
	//normal s;
	//cout << s.mean() << endl;
	//	int n = 1 << order;
	//	//Hilbert_curve curve(n);
	//	cout << n;
	//	Zig_zag_curve curve(pow(2, order),3);
	//	curve.get_points_for_curve();
	//	std::cout << "Hilbert curve, order=" << order << '\n';
	//	/*std::vector<std::string> lines = curve.draw_curve();
	//	for (auto &line : lines) {
	//		std::cout << line.c_str() << '\n';
	//	}*/
	//	for (auto& point : curve.get_points())
	//		std::cout << '(' << point.x << ',' << point.y << ')' << '\n';
	//	std::cout << '\n';
	//}
	//initial_prob_img *img = new initial_prob_img(32, NORMAL, 1, 2);
	//initial_prob_img *img = new initial_prob_img(32, RAYLEIGH, 1, 2);
	/*initial_prob_img *img = new initial_prob_img("D:\\generated_image.txt", 32, RAYLEIGH, 1, 2);*/
	//initial_prob_img *img = new initial_prob_img("D:\\generated_image.txt", 32, NORMAL, 1, 2);
	initial_prob_img *img = new initial_prob_img("C:\\Users\\anastasya\\Desktop\\data_FR.txt");
	//initial_prob_img *img = new initial_prob_img("C:\\Users\\anastasya\\Desktop\\data_test.txt");
	
	shared_ptr<initial_prob_img> ptr_img(img);
	//mixture_handler t(ptr_img->get_m_image(), 5, 0.001);
	mixture_handler t(ptr_img, 5, 0.001);
	//t.mixture_optimal_redraw_opMP_V2();
	//t.kolmogorov_optimal_redraw_opMP();
	t.q_tree_optimal_redraw_opMP();
	t.printInformation();
	t.draw_graphics();
	//quad_tree_handler bbbb = quad_tree_handler(ptr_img, 32);
	//bbbb.set_probabilities(0, 0);
	//bbbb.bottom_up_pass();
	//bbbb.up_down_pass();
	////bbbb.split_image_by_summ();
	//bbbb.split_image_by_vote();
	//bbbb.create_splitted_img();
	//bbbb.draw_graphics();
	//path p{ "C:\\Windows\\System" };
	//std::cout << p.string() << '\n';
	//mix_img_obj("C:\\Users\\anastasya\\Desktop\\data_test.txt");
	return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
