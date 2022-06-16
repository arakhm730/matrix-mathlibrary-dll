// MathLibrary.h - Contains declarations of math functions
#pragma once


#ifdef MATHLIBRARY_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif


#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <iterator>
#include <sstream>

#  pragma warning( push )
#  pragma warning( disable: 4251 )

/*class Matrix;
class DiagonalMatrix;
class UpperTriangularMatrix;
class LowerTriangularMatrix;
class SymmetricalMatrix;
*/



class MATHLIBRARY_API Matrix {
private:
	int Row, Col;
	std::vector <std::vector <double>> mtr;



	double getDeterminant(const std::vector<std::vector<double>> vect);
	std::vector<std::vector<double>> getTranspose(const std::vector<std::vector<double>> matrix1);
	std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>> vect);
	std::vector<std::vector<double>> getInverse(const std::vector<std::vector<double>> vect);
	int compute_rank(std::vector<std::vector<double>> A);

public:
	Matrix();
	Matrix(std::vector <std::vector <double>> data);
	Matrix(const int&, const int&);
	virtual ~Matrix();
	Matrix(const Matrix&);
	std::vector<std::vector<double>> GetMtrData() const;
	int GetRow() const;
	int GetCol() const;

	void SetRow(const int&);
	void SetCol(const int&);
	void SetMtr(std::vector<std::vector<double>>);

	// Лабораторная работа №2

	double determinant();
	double trace() const;
	double vecNormMax() const;
	double vecNormEuclid() const;
	double matrixNorm() const;

	Matrix transpose();

	// Лабораторная работа №5
	Matrix center();
	Matrix scaling();
	std::vector<Matrix> pca(unsigned int pc);

	int rank();
	Matrix inverse();

	std::ifstream& load_from(std::ifstream&);
	std::ofstream& save_to(std::ofstream&);

	friend class IdentityMatrix;
	friend class DiagonalMatrix;
	friend class UpperTriangularMatrix;
	friend class LowerTriangularMatrix;
	friend class SymmetricalMatrix;

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);

};




class MATHLIBRARY_API IdentityMatrix : public Matrix {
public:

	IdentityMatrix(const int& dim) : Matrix::Matrix(dim, dim) {
		for (int i = 0; i < Row; i++) for (int j = 0; j < Col; j++)
			i == j ? this->mtr.at(i).at(j) = 1 : this->mtr.at(i).at(j) = 0;
	}
	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);
};


class MATHLIBRARY_API DiagonalMatrix : public Matrix {
public:

	DiagonalMatrix(const Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);
};


class MATHLIBRARY_API UpperTriangularMatrix : public Matrix {
public:

	UpperTriangularMatrix(const Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);
};


class MATHLIBRARY_API LowerTriangularMatrix : public Matrix {
public:

	LowerTriangularMatrix(const Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);
};


class MATHLIBRARY_API SymmetricalMatrix : public Matrix {
public:
	SymmetricalMatrix(const Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix AdamarMult(Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator+ (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator- (Matrix&, Matrix&);

	friend MATHLIBRARY_API Matrix operator* (const double&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, const double&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, Matrix&);
	friend MATHLIBRARY_API Matrix operator* (Matrix&, std::vector<double>&);
	friend MATHLIBRARY_API Matrix operator* (std::vector<double>&, std::vector<double>&);
	//friend std::ifstream& operator>> (std::ifstream& ifs, Matrix<double>& fill);
	friend MATHLIBRARY_API std::ofstream& operator<< (std::ofstream&, const Matrix&);
	friend MATHLIBRARY_API std::ostream& operator<< (std::ostream&, const Matrix&);
	friend MATHLIBRARY_API std::istream& operator>> (std::istream&, Matrix&);

	friend MATHLIBRARY_API double swing(Matrix&, unsigned int&);
	friend MATHLIBRARY_API double leftover(Matrix&, unsigned int&);
	friend MATHLIBRARY_API std::vector<double> dispersions(Matrix&, Matrix&);
};

extern "C" MATHLIBRARY_API std::ifstream& operator>>(std::ifstream&, Matrix&);