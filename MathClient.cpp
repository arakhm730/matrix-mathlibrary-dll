// MathClient.cpp : Client app for MathLibrary DLL.
//#include "pch.h"
#include "MathLibrary.h"

int main() {
	setlocale(LC_ALL, "RUS");

	std::ifstream ifs("data.txt");
	Matrix txt;
	ifs >> txt;
	ifs.close();

	std::ofstream ofs("out.txt");
	ofs.clear();
	ofs << txt;
	ofs.close();



	std::ofstream out("numbers", std::ios::out | std::ios::binary);
	out.clear();


	txt.save_to(out);
	out.close();


	Matrix ot({ {1, 2, 3, 4}, {-5, -1, 2, 10}, {-1, 7, 1.5, 2}, {0, 1, 1, 0} });
	//std::vector<Matrix> pc = ot.pca(3);
	//for (auto& elm : pc)	std::cout << elm;

	std::vector<double> a = { 1, 2, 3, 4 };
	std::cout << ot << std::endl;

	//std::cout << "after multiplying:\n";
	//ot = ot * a;
	//std::cout << ot << std::endl;
	std::cout << "Matrix ot after centering:\n";
	ot.center();
	std::cout << ot << std::endl;

	std::cout << "Matrix ot after scaling:\n";
	ot.scaling();
	std::cout << ot << std::endl;






	std::cout << "Matrix txt after centering:\n";
	txt.center();
	std::cout << txt << std::endl;

	std::cout << "Matrix txt after scaling (Pre-PCA Look):\n";
	txt.scaling();
	std::cout << txt << std::endl;
	std::vector<Matrix> pc = txt.pca(5);

	std::cout << "\n\n\n////// PCA Method: 1st is matrix P, 2nd is matrix T, 3rd is matrix E //////\n";
	for (auto& elm : pc)	std::cout << elm;
	std::cout << "\n\n\n";

	Matrix Pt = pc[0].transpose();
	Matrix tp = pc[1] * Pt;
	Matrix x_like = tp + pc[2];

	std::cout << "x_like matrix (T*P' + E):\n" << x_like << std::endl;

	unsigned int row = 4;
	std::cout << "Swing of row #" << row << ":\n" << swing(pc[1], row) << std::endl;
	std::cout << "Dispersions TRV and ERV respectively:\n";
	std::vector<double> disp = dispersions(txt, pc[2]);
	for (auto& elm : disp)	std::cout << elm << " ";
	std::cout << std::endl;
	std::ifstream in("numbers", std::ios::in | std::ios::binary);
	Matrix ii(4, 4);
	ii.load_from(in);
	std::cout << ii << std::endl;

	in.close();
	std::cout << "slatt?" << std::endl;

	return 0;
}