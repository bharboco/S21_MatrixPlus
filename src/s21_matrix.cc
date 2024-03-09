#include "s21_matrix_oop.h"

// Методы

S21Matrix::S21Matrix() : rows_(1), cols_(1) {
  AlocMatrix(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  AlocMatrix(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  matrix_ = other.matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    S21Matrix::DelMatrix(matrix_);
  }
}

void S21Matrix::DelMatrix(double** matrix) {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

void S21Matrix::AlocMatrix(double*** matrix, int rows, int cols) {
  if (rows < 1)
    throw std::invalid_argument("неккоректное значение строк -> " +
                                std::to_string(rows));
  if (cols < 1)
    throw std::invalid_argument("неккоректное значение столбцов -> " +
                                std::to_string(cols));
  *matrix = new double*[rows]();
  for (int i = 0; i < rows; i++) (*matrix)[i] = new double[cols]();
}

// Операции

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7)
        return false;  // 0.0000001
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Различная размерность матриц");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Различная размерность матриц");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double& num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = (*this)(i, j) * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_)
    throw std::invalid_argument(
        "Неккоректная матрица : число столбцов первой матрицы не равно числу "
        "строк второй матрицы.");
  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        res.matrix_[i][j] += (*this)(i, k) * other.matrix_[k][j];
      }
    }
  }
  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix trans(cols_, rows_);
  for (int i = 0; i < trans.rows_; i++) {
    for (int j = 0; j < trans.cols_; j++) {
      trans.matrix_[i][j] = matrix_[j][i];
    }
  }
  return trans;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::invalid_argument("Матрица не является квадратной");
  S21Matrix calc(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      calc.matrix_[i][j] = pow(-1, (i + j)) * Minor(i, j).Determinant();
    }
  }
  return calc;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_)
    throw std::invalid_argument("Матрица не является квадратной");
  double rezult = 0;
  if (rows_ == 1) {
    rezult = (*this)(0, 0);
  } else if (rows_ == 2) {
    rezult = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    int sign = 1;
    for (int i = 0; i < cols_; i++) {
      rezult += (*this)(0, i) * Minor(0, i).Determinant() * sign;
      sign = -sign;
    }
  }
  return rezult;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_)
    throw std::invalid_argument("Матрица не является квадратной");
  double d = Determinant();
  if (fabs(d) < 1e-7) throw std::invalid_argument("Матрица равна 0");
  S21Matrix inverse(rows_, cols_);
  inverse = CalcComplements().Transpose();
  inverse.MulNumber(1 / d);
  return inverse;
}

// мутатор и ацессор
int S21Matrix::getrows() const { return rows_; }

int S21Matrix::getcols() const { return cols_; }

void S21Matrix::setrows(int rows) {
  if (rows == rows_) return;
  double** new_matrix;
  AlocMatrix(&new_matrix, rows, cols_);
  for (int i = 0; i < rows && i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      new_matrix[i][j] = matrix_[i][j];
    }
  }
  if (rows > rows_) {
    for (int i = rows_; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        new_matrix[i][j] = 0;
      }
    }
  }
  DelMatrix(matrix_);
  rows_ = rows;
  matrix_ = new_matrix;
}

void S21Matrix::setcols(int cols) {
  if (cols == cols_) return;
  double** new_matrix;
  AlocMatrix(&new_matrix, rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols && j < cols_; j++) {
      new_matrix[i][j] = matrix_[i][j];
    }
  }
  if (cols > cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = cols_; j < cols; j++) {
        new_matrix[i][j] = 0;
      }
    }
  }
  DelMatrix(matrix_);
  cols_ = cols;
  matrix_ = new_matrix;
}

// оператор
S21Matrix S21Matrix::operator+(S21Matrix& other) {
  S21Matrix rezult = *this;
  rezult.SumMatrix(other);
  return rezult;
}
S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix rezult = *this;
  rezult.SubMatrix(other);
  return rezult;
}
S21Matrix S21Matrix::operator*(S21Matrix& other) {
  S21Matrix rezult = *this;
  rezult.MulMatrix(other);
  return rezult;
}
S21Matrix S21Matrix::operator*(const double& num) {
  S21Matrix rezult = *this;
  MulNumber(num);
  return rezult;
}
S21Matrix S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix(other.rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}
bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}
void S21Matrix::operator+=(const S21Matrix& other) { return SumMatrix(other); }
void S21Matrix::operator-=(const S21Matrix& other) { return SubMatrix(other); }
void S21Matrix::operator*=(const S21Matrix& other) { return MulMatrix(other); }
void S21Matrix::operator*=(const double num) { return MulNumber(num); }
double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || j < 0) {
    throw std::invalid_argument("Недействительные значения строк и столбцов");
  }
  return matrix_[i][j];
}

// tools
S21Matrix S21Matrix::Minor(int i_rows, int j_cols) {
  if (i_rows < 0 || i_rows >= rows_ || j_cols < 0 || j_cols >= cols_) {
    throw std::invalid_argument("Недействительные значения строк и столбцов");
  }
  int rows = rows_ - 1, cols = cols_ - 1;
  if (rows == 0) rows = 1;
  if (cols == 0) cols = 1;
  S21Matrix minor(rows, cols);
  int minor_rows = 0;
  for (int i = 0; i < rows_; i++) {
    int minor_cols = 0;
    if (i == i_rows) continue;
    for (int j = 0; j < cols_; j++) {
      if (j == j_cols) continue;
      minor.matrix_[minor_rows][minor_cols] = (*this)(i, j);
      minor_cols++;
    }
    minor_rows++;
  }
  return minor;
}