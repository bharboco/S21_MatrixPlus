#ifndef CC1_S21_MATRIXPLUS_H
#define CC1_S21_MATRIXPLUS_H

#include <cmath>
#include <iostream>

class S21Matrix {
 public:
  // Методы
  S21Matrix();  // конструктор по заранее заданной размерности
  S21Matrix(int rows, int cols);  //конструктор по параметрам
  S21Matrix(const S21Matrix& other);  //  копирование
  S21Matrix(S21Matrix&& other);       // перенос
  ~S21Matrix();                       // деконструктор

  void AlocMatrix(double*** matrix, int rows, int cols);
  void DelMatrix(double** matrix_);

  //Операции
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double& num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // tools
  S21Matrix Minor(int i_rows, int j_cols);
  S21Matrix Check(int i, int j);

  // accessor и mutator
  int getrows() const;
  int getcols() const;
  void setrows(int rows);
  void setcols(int cols);

  //Оператор
  S21Matrix operator+(S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(S21Matrix& other);
  S21Matrix operator*(const double& num);
  S21Matrix operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other) const;
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double num);
  double& operator()(int i, int j);

 private:
  int rows_, cols_;  //строки , столбцы
  double** matrix_;
};

#endif