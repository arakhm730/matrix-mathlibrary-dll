// MathLibrary.cpp : Defines the exported functions for the DLL.
#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier
#include <utility>
#include <limits.h>
#include "MathLibrary.h"



const double eps = 10e-8;

Matrix::Matrix(const int& rows, const int& cols) : Row(rows), Col(cols) {
	try
	{
		if (rows < 0 || cols < 0) throw(123);
		mtr.assign(rows, std::vector<double>(cols));
	}
	catch (const int) {
		std::cout << "Неправильно задана размерность матрицы!!" << std::endl;
	}
}

Matrix::Matrix() : mtr({ {} }), Row(0), Col(0) {}

Matrix::Matrix(std::vector <std::vector <double>> data) : mtr(data), Row(data.size()), Col(data.at(0).size()) {}

Matrix::Matrix(const Matrix& sample) {
	int rows = sample.GetRow();
	int cols = sample.GetCol();

	try
	{
		if (rows < 0 || cols < 0) throw('e');
		mtr.assign(rows, std::vector<double>(cols));
	}
	catch (const char) {
		std::cout << "Неправильно задана размерность матрицы!" << std::endl;
	}

	this->Row = rows;
	this->Col = cols;



	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			mtr.at(i).at(j) = sample.mtr.at(i).at(j);
}


Matrix::~Matrix() = default;

int Matrix::GetRow() const {
	return Row;
}

int Matrix::GetCol() const {
	return Col;
}


void Matrix::SetRow(const int& row)
{
	Row = row;
}

void Matrix::SetCol(const int& col)
{
	Col = col;
}


void Matrix::SetMtr(std::vector<std::vector<double>> data)
{
	mtr = data;
}


std::vector<std::vector<double>> Matrix::GetMtrData() const {
	return mtr;
}


double Matrix::trace() const {
	try {
		if (this->GetCol() != this->GetRow()) throw (2.34);
		else {
			double tr = 0.0;
			for (int i = 0; i < this->GetRow(); i++) tr += this->mtr.at(i).at(i);
			return tr;
		}
	}
	catch (double) { std::cout << "Cannot calculate trace of non-square matrix!" << std::endl; }
}


double Matrix::vecNormMax() const {
	try {
		if (this->Row != 1) throw(-5);
		else {
			double max = 0.0;
			for (int j = 0; j < this->Col; j++)
				if ((double)this->mtr.at(0).at(j) >= max) max = this->mtr.at(0).at(j);
			return max;
		}
	}
	catch (int) { std::cout << "Cannot calculate vector norm of non-vector object!" << std::endl; }
}


double Matrix::vecNormEuclid() const {
	try {
		if (this->Row != 1) throw(-5);
		else {
			double squares = 0.0;
			for (int j = 0; j < this->Col; j++) squares += pow(this->mtr.at(0).at(j), 2);
			return (double)pow(squares, 0.5);
		}
	}
	catch (int) { std::cout << "Cannot calculate vector norm of non-vector object!" << std::endl; }
}


double scalarProduct(const Matrix& first, const Matrix& second) {
	try {
		if (((first.GetRow() == second.GetRow() && second.GetRow() == 1) && first.GetCol() == second.GetCol())) {
			double product = 0.0;
			for (int i = 0; i < first.GetCol(); i++) product += first.GetMtrData().at(0).at(i) * second.GetMtrData().at(0).at(i);
			return product;
		}
		else { throw (1); }
	}
	catch (int) { std::cout << "Vectors have inequal sizes, cannot calculate scalar product!" << std::endl; }
}


double Matrix::matrixNorm() const {
	double squares = 0;
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			squares += pow(this->mtr.at(i).at(j), 2);
	return pow((double)squares, 0.5);
}

double vectors_angle(const Matrix& first, const Matrix& second) {
	double scalar = scalarProduct(first, second);
	return acos(scalar / (first.vecNormEuclid() * second.vecNormEuclid()));
}

double det(const Matrix& inst, const int n) {
	try {
		if (inst.GetCol() != inst.GetRow()) throw (2.34);
		else {
			double deter = 0;
			Matrix submatrix(n, n);
			if (n == 2)
				return ((inst.GetMtrData().at(0).at(0) * inst.GetMtrData().at(1).at(1)) - (inst.GetMtrData().at(1).at(0) * inst.GetMtrData().at(0).at(1)));
			else {
				for (int x = 0; x < n; x++) {
					int subi = 0;
					for (int i = 1; i < n; i++) {
						int subj = 0;
						for (int j = 0; j < n; j++) {
							if (j == x)
								continue;
							submatrix.GetMtrData().at(subi).at(subj) = inst.GetMtrData().at(i).at(j);
							subj++;
						}
						subi++;
					}
					deter = deter + (pow(-1, x) * inst.GetMtrData().at(0).at(x) * det(submatrix, n - 1));
				}
			}
			return deter;
		}
	}
	catch (double) { std::cout << "Cannot calculate determinant of non-square matrix!" << std::endl; }
}

double Matrix::determinant() {
	Matrix copy(this->Row, this->Col);
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			copy.mtr[i][j] = this->mtr[i][j];
	return getDeterminant(copy.mtr);
}



Matrix Matrix::transpose() {
	std::vector <std::vector<double>> copy(this->Col, std::vector<double>());
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++) {
			copy.at(j).push_back(this->mtr.at(i).at(j));
		}
	Matrix c(copy);
	return c;
}


double Matrix::getDeterminant(const std::vector<std::vector<double>> vect) {
	if (vect.size() != vect[0].size()) {
		throw std::runtime_error("Matrix is not quadratic");
	}
	int dimension = vect.size();

	if (dimension == 0) {
		return 1;
	}

	if (dimension == 1) {
		return vect[0][0];
	}

	//Formula for 2x2-matrix
	if (dimension == 2) {
		return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
	}

	double result = 0;
	int sign = 1;
	for (int i = 0; i < dimension; i++) {

		//Submatrix
		std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double>(dimension - 1));
		for (int m = 1; m < dimension; m++) {
			int z = 0;
			for (int n = 0; n < dimension; n++) {
				if (n != i) {
					subVect[m - 1][z] = vect[m][n];
					z++;
				}
			}
		}

		//recursive call
		result = result + sign * vect[0][i] * getDeterminant(subVect);
		sign = -sign;
	}

	return result;
}


std::vector<std::vector<double>> Matrix::getTranspose(const std::vector<std::vector<double>> matrix1) {

	//Transpose-matrix: height = width(matrix), width = height(matrix)
	std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double>(matrix1.size()));

	//Filling solution-matrix
	for (size_t i = 0; i < matrix1.size(); i++) {
		for (size_t j = 0; j < matrix1[0].size(); j++) {
			solution[j][i] = matrix1[i][j];
		}
	}
	return solution;
}

std::vector<std::vector<double>> Matrix::getCofactor(const std::vector<std::vector<double>> vect) {
	if (vect.size() != vect[0].size()) {
		throw std::runtime_error("Matrix is not quadratic");
	}

	std::vector<std::vector<double>> solution(vect.size(), std::vector<double>(vect.size()));
	std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double>(vect.size() - 1));

	for (std::size_t i = 0; i < vect.size(); i++) {
		for (std::size_t j = 0; j < vect[0].size(); j++) {

			int p = 0;
			for (size_t x = 0; x < vect.size(); x++) {
				if (x == i) {
					continue;
				}
				int q = 0;

				for (size_t y = 0; y < vect.size(); y++) {
					if (y == j) {
						continue;
					}

					subVect[p][q] = vect[x][y];
					q++;
				}
				p++;
			}
			solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
		}
	}
	return solution;
}


std::vector<std::vector<double>> Matrix::getInverse(const std::vector<std::vector<double>> vect) {
	if (getDeterminant(vect) == 0) {
		throw std::runtime_error("Determinant is 0");
	}

	double d = 1.0 / getDeterminant(vect);
	std::vector<std::vector<double>> solution(vect.size(), std::vector<double>(vect.size()));

	for (size_t i = 0; i < vect.size(); i++) {
		for (size_t j = 0; j < vect.size(); j++) {
			solution[i][j] = vect[i][j];
		}
	}

	solution = getTranspose(getCofactor(solution));

	for (size_t i = 0; i < vect.size(); i++) {
		for (size_t j = 0; j < vect.size(); j++) {
			solution[i][j] *= d;
		}
	}

	return solution;
}


Matrix Matrix::inverse() {
	Matrix copy(this->mtr);
	Matrix inversed(getInverse(copy.GetMtrData()));
	return inversed;
}

const double EPS = 1E-9;


int Matrix::compute_rank(std::vector<std::vector<double>> A) {
	int n = A.size();
	int m = A[0].size();

	int rg = 0;
	std::vector<bool> row_selected(n, false);
	for (int i = 0; i < m; ++i) {
		int j;
		for (j = 0; j < n; ++j) {
			if (!row_selected[j] && abs(A[j][i]) > EPS)
				break;
		}

		if (j != n) {
			++rg;
			row_selected[j] = true;
			for (int p = i + 1; p < m; ++p)
				A[j][p] /= A[j][i];
			for (int k = 0; k < n; ++k) {
				if (k != j && abs(A[k][i]) > EPS) {
					for (int p = i + 1; p < m; ++p)
						A[k][p] -= A[j][p] * A[k][i];
				}
			}
		}
	}
	return rg;
}


Matrix Matrix::center() {
	for (int j = 0; j < Col; ++j) {
		double m = 0.0;
		for (int i = 0; i < Row; ++i)	m += mtr[i][j];
		for (int i = 0; i < Row; ++i)	mtr[i][j] -= (double)m / Row;
	}
	return (*this);
}


Matrix Matrix::scaling() {
	for (int j = 0; j < Col; ++j) {
		double m = 0.0;
		for (int i = 0; i < Row; ++i)	m += mtr[i][j];
		for (int i = 0; i < Row; ++i)	mtr[i][j] -= (double)m / Row;

		double square_sum = 0;
		for (int i = 0; i < Row; ++i)	square_sum += mtr[i][j] * mtr[i][j];
		if (square_sum == 0)
			std::cout << "standard deviation is zero, no scaling executed in column #" << j << std::endl;
		else
			for (int i = 0; i < Row; ++i)	mtr[i][j] /= sqrt(1.0 / (Row - 1) * square_sum);
	}
	return (*this);
}




int Matrix::rank() {
	return compute_rank(this->mtr);
}


Matrix AdamarMult(Matrix& first, Matrix& second) {
	try {
		if (first.Row != second.Row || first.Col != second.Col) throw("help");

		Matrix c(first.Row, first.Col);
		for (int i = 0; i < first.Row; i++)
			for (int j = 0; j < first.Col; j++)
				c.mtr.at(i).at(j) = first.mtr.at(i).at(j) * second.mtr.at(i).at(j);
		return c;

	}
	catch (const char*) {
		std::cout << "Incompatible matrix sizes!" << std::endl;
	}
}


Matrix operator+ (Matrix& first, Matrix& second) {
	if (first.Row == second.Row && first.Col == second.Col) {
		Matrix c(first.Row, first.Col);
		for (int i = 0; i < first.Row; i++)
			for (int j = 0; j < first.Col; j++)
				c.mtr.at(i).at(j) = first.mtr.at(i).at(j) + second.mtr.at(i).at(j);
		return c;
	}
	else std::cout << "Incompatible matrix sizes!" << std::endl;

}


Matrix operator- (Matrix& first, Matrix& second) {
	if (first.Row == second.Row && first.Col == second.Col) {
		Matrix c(first.Row, first.Col);
		for (int i = 0; i < first.Row; i++)
			for (int j = 0; j < first.Col; j++)
				c.mtr.at(i).at(j) = first.mtr.at(i).at(j) - second.mtr.at(i).at(j);
		return c;
	}
	else std::cout << "Incompatible matrix sizes!" << std::endl;

}


Matrix operator*(const double& arg, Matrix& inst) {
	Matrix multi(inst.GetRow(), inst.GetCol());
	for (int i = 0; i < inst.GetRow(); i++)
		for (int j = 0; j < inst.GetCol(); j++)
			multi.mtr.at(i).at(j) = inst.mtr.at(i).at(j) * arg;
	return multi;
}


Matrix operator*(Matrix& inst, const double& arg) {
	Matrix multi(inst.GetRow(), inst.GetCol());
	for (int i = 0; i < inst.GetRow(); i++)
		for (int j = 0; j < inst.GetCol(); j++)
			multi.mtr.at(i).at(j) = inst.mtr.at(i).at(j) * arg;
	return multi;
}


Matrix operator*(Matrix& inst1, Matrix& inst2) {
	if (inst1.GetCol() == inst2.GetRow()) {

		Matrix mult(inst1.GetRow(), inst2.GetCol());
		for (int i = 0; i < inst1.GetRow(); i++) {
			for (int j = 0; j < inst2.GetCol(); j++) {
				double elm = 0.0;
				for (int k = 0; k < inst1.GetCol(); k++)
					elm += inst1.mtr.at(i).at(k) * inst2.mtr.at(k).at(j);

				mult.mtr.at(i).at(j) = elm;
			}
		}
		return mult;
	}
	else std::cout << "Incompatible matrices size!" << std::endl;

}


Matrix operator* (std::vector<double>& vec, Matrix& matrix) {
	std::vector<std::vector<double>> mat;
	mat.emplace_back(vec);
	Matrix mt = Matrix(mat);
	return (mt * matrix);
}


Matrix operator* (Matrix& matrix, std::vector<double>& vec) {
	std::vector<std::vector<double>> mat;
	mat.emplace_back(vec);
	Matrix _matrix = Matrix(mat);
	_matrix = _matrix.transpose();
	return (matrix * _matrix);
}

Matrix operator* (std::vector<double>& first, std::vector<double>& second) {
	std::vector<std::vector<double>> vecs1;
	vecs1.emplace_back(first);
	std::vector<std::vector<double>> vecs2;
	vecs2.emplace_back(second);
	Matrix m1 = Matrix(vecs1);
	m1 = m1.transpose();
	Matrix m2 = Matrix(vecs2);
	return (m1 * m2);
}


double vecNorm(std::vector<double> vec) {
	double sum = 0;
	for (auto& elm : vec)	sum += (double)pow(elm, 2);
	return sqrt(sum);
}


std::ostream& operator<< (std::ostream& out, const Matrix& instance) {
	out << "Current matrix data:\n";
	int Row = instance.GetRow();
	int Col = instance.GetCol();
	for (int i = 0; i < Row; i++)
	{
		if (i == 0) out << "/ ";
		else if (i == Row - 1) out << "\\ ";
		else out << "| ";
		for (int j = 0; j < Col; j++) {
			out << round(instance.mtr.at(i).at(j) * 1000) / 1000 << std::setw(1) << '\t';
			//if (j != Col - 1) out
		}


		if (i == 0) out << " \\";
		else if (i == Row - 1) out << " /";
		else out << " |";
		out << "\n";
	}
	return out;
}


std::istream& operator>> (std::istream& in, Matrix& dest) {
	int Row = dest.GetRow();
	int Col = dest.GetCol();
	std::cout << "Пожалуйста, введите " << Row << " строк(-и), в каждой из которых по " << Col << " элемента(-ов), разделяя элементы в строке пробелами, а строки отделяя новой строкой:\n";
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			in >> dest.mtr.at(i).at(j);
	return in;
}


std::ofstream& operator<<(std::ofstream& ofs, const Matrix& data) {
	for (int i = 0; i < data.GetRow(); i++) {
		for (int j = 0; j < data.GetCol(); j++) {
			ofs << data.GetMtrData().at(i).at(j);
			j != data.GetCol() - 1 ? (ofs << '\t') : (ofs << '\n');
		}
	}
	return ofs;
}




DiagonalMatrix::DiagonalMatrix(const Matrix& inst) : Matrix::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i == j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}


UpperTriangularMatrix::UpperTriangularMatrix(const Matrix& inst) : Matrix::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i >= j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}


LowerTriangularMatrix::LowerTriangularMatrix(const Matrix& inst) : Matrix::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i <= j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}


SymmetricalMatrix::SymmetricalMatrix(const Matrix& inst) : Matrix::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i > j ? this->mtr.at(i).at(j) = inst.mtr.at(j).at(i) : this->mtr.at(i).at(j) = inst.mtr.at(i).at(j);
}


std::ifstream& Matrix::load_from(std::ifstream& file) {
	if (!file) {
		std::cout << "Cannot open file.";
		return file;
	}
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			file.read(reinterpret_cast<char*>(&mtr[i][j]), sizeof(double));

	return file;
}


std::ofstream& Matrix::save_to(std::ofstream& bin) {
	if (!bin) {
		std::cout << "Cannot open file.";
		return bin;
	}

	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			bin.write(reinterpret_cast<char*>(&mtr[i][j]), sizeof(double));

	return bin;

}

#define minimum(a, b) ((a) < (b) ? a : b)
#define maximum(a, b) ((a) > (b) ? a : b)
std::vector<Matrix> Matrix::pca(unsigned int pc) {
	try
	{
		if (pc < 1 || pc > minimum(this->Row, this->Col))	throw(std::runtime_error("Not appropriate number for pca algorithm!"));
		Matrix E = (*this);
		Matrix P(std::vector<std::vector<double>>(1, std::vector<double>(minimum(Row, Col), 0)));
		Matrix T_c(std::vector<std::vector<double>>(1, std::vector<double>(maximum(Row, Col), 0)));
		Matrix E_t = E.transpose();
		for (unsigned int h = 0; h < pc; ++h) {

			std::vector<double> t = E_t.mtr.at(h);

			std::vector<double> p;
			std::vector<double> t_old;
			std::vector<double> d;
			do {
				Matrix tm0 = t * E;
				Matrix tm = tm0 * (1 / pow(vecNorm(t), 2));
				p = tm.mtr.at(0);
				for (auto& elm : p)	elm /= vecNorm(p);

				t_old = t;

				Matrix t0 = E * p;
				t0 = t0 * (1 / pow(vecNorm(p), 2));
				t0 = t0.transpose();
				t = t0.mtr.at(0);
				d.clear();
				for (size_t i = 0; i < t.size(); ++i)	d.emplace_back(t_old.at(i) - t.at(i));
			} while (vecNorm(d) > eps);
			Matrix e0 = t * p;
			E = E - e0;
			P.mtr.emplace_back(p);
			P.Row++;
			T_c.mtr.emplace_back(t);
			T_c.Row++;
		}
		P = Matrix(P.mtr);
		P = P.transpose();
		for (auto& elm : P.mtr)	elm.erase(elm.begin()); // Удалить нули в начале
		P = Matrix(P.mtr);

		T_c = Matrix(T_c.mtr);
		T_c = T_c.transpose();
		for (auto& elm : T_c.mtr)	elm.erase(elm.begin()); // Удалить нули в начале
		T_c = Matrix(T_c.mtr);

		std::vector<Matrix> ret;
		ret.emplace_back(P);
		ret.emplace_back(T_c);
		ret.emplace_back(E);
		Matrix P_t = P.transpose();
		Matrix I = P_t * P;
		std::cout << "\nP_t * P:\n" << I << "\n\n\n";
		return ret;
	}
	catch (std::runtime_error& f) { std::cout << f.what() << std::endl;	return std::vector<Matrix>(); }
}

double swing(Matrix& T_scores, unsigned int& row) {
	Matrix T_t = T_scores.transpose();
	Matrix Multi = T_t * T_scores;

	// Inversing diagonal matrix:
	for (int i = 0; i < Multi.Row; ++i) {
		for (int j = 0; j < Multi.Col; ++j) {
			if (j != i)	Multi.mtr[i][j] = 0;
			else
			{
				Multi.mtr[i][j] = 1 / Multi.mtr[i][j];
			}
		}
	}

	Matrix pre_answer = T_scores.mtr.at(row) * Multi;
	Matrix answer = pre_answer * T_scores.mtr.at(row);
	return answer.mtr.at(0).at(0);
}

double leftover(Matrix& E, unsigned int& row) {
	double answer = 0.0;
	for (int j = 0; j < E.Col; ++j)	answer += pow(E.mtr[row][j], 2);
	return answer;
}

std::vector<double> dispersions(Matrix& X, Matrix& E) {
	double lo_sum = 0;
	for (unsigned int i = 0; i < E.Row; ++i)	lo_sum += leftover(E, i);
	lo_sum /= E.Row;

	double TRV = lo_sum / E.Col;

	double e = 0.0;
	double x = 0.0;

	for (int i = 0; i < X.Row; ++i) {
		for (int j = 0; j < X.Col; ++j)	x += pow(X.mtr[i][j], 2);
	}
	x /= X.Col;
	x /= X.Row;

	double ERV = 1 - TRV / x;

	std::vector<double> answer;
	answer.push_back(TRV);
	answer.push_back(ERV);
	return answer;

}

std::ifstream& operator>>(std::ifstream& ifs, Matrix& fill) {
	try {
		if (ifs.is_open()) {
			std::vector<std::vector<double>> numbers;

			std::string temp;
			int size = 0; int size1 = 0;
			int k = 0;
			std::cout << "Imported file contents:\n----------------------------------------------------\n";
			while (std::getline(ifs, temp)) {
				std::istringstream buffer(temp);
				std::vector<double> line;
				for (std::string s; buffer >> s;) line.push_back(std::stod(s));
				size = line.size();
				if (size != size1 && k != 0) throw(std::runtime_error("Cannot build matrix!"));
				size1 = line.size();
				for (auto elm : line) std::cout << elm << " ";
				std::cout << std::endl;

				numbers.push_back(line);
				k++;
			}
			fill.SetRow(numbers.size());
			fill.SetCol(numbers[0].size());
			fill.SetMtr(numbers);
			std::cout << "----------------------------------------------------\n";
		}

		else { std::cout << "File wasn't opened!" << std::endl; }
	}
	catch (std::runtime_error&) { std::cout << "\nFix the text file contents and make it look like 2D matrix.\n"; };
	return ifs;
}