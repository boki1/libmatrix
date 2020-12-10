#include <algorithm>
#include <iostream>
#include <vector>

namespace matrix {

class MatrixException : public std::exception {};

class UndefinedBahaviourMatrixException : public MatrixException {
private:
  const std::string context;

public:
  UndefinedBahaviourMatrixException(std::string _context) : context(_context) {}
  const std::string &Why() { return this->context; }
};

template <typename T> class Matrix {
public:
  Matrix();

  explicit Matrix(int d, int x = Matrix::N, int y = Matrix::N);

  Matrix(const Matrix &other)
      : data(other.data), dimensions(other.dimensions) {}

  Matrix(const Matrix &&other)
      : data(other.data), dimensions(other.dimensions) {
    data.clear();
    set_dimensions(0, 0);
  }

  Matrix(std::vector<std::vector<T>> data);

  ~Matrix() = default;

  Matrix &normalize_row(int row, double coeff);

  T *get_row(int i) const;

  Matrix &swap_row(int i, int j);

  Matrix &save_row(int i, T *values);

  Matrix &multadd_rows(int i, int j, T val);

  Matrix<T> reduce() const;

  void set_empty() {
    for (int i = 0; i < Matrix::N * Matrix::N; ++i) {
      data.push_back(0);
    }
    dimensions = {0, 0};
  }

  void set_diagonal(int v) {
    for (int i = 0; i < y_dim(); ++i) {
      for (int j = 0; j < x_dim(); ++j) {
        if (i == j)
          this->store_at(j, i, v);
      }
    }
  }

  int same_dimensions_as(const Matrix &other) {
    if (this->dimensions.first == other.dimensions.first &&
        this->dimensions.second == other.dimensions.second) {
      return 1;
    }
    return 0;
  }

  Matrix &transposed() const {
    Matrix *t = new Matrix;
    t->set_dimensions(this->y_dim(), this->x_dim());
    int val;
    for (int y = 0; y < this->dimensions.second; ++y) {
      for (int x = 0; x < this->dimensions.first; ++x) {
        val = this->at(x, y);
        t->store_at(y, x, val);
      }
    }
    return *t;
  }

  int index(int x, int y) const {
    int idx = y * dimensions.first + x;
    return idx;
  }

  T at(int x, int y) const {
    int idx = index(x, y);
    return this->data[idx];
  }

  void store_at(int x, int y, T val) {
    int idx = index(x, y);
    this->data[idx] = val;
  }

  int x_dim() const { return this->dimensions.first; }

  int y_dim() const { return this->dimensions.second; }

  void set_dimensions(int x, int y) { this->dimensions = {x, y}; }

  std::pair<int, int> get_dimensions() const { return this->dimensions; }

  std::vector<T> get_data() const { return this->data; }

  static const unsigned N = 50;

  template <typename U>
  friend int operator==(matrix::Matrix<U> &This, matrix::Matrix<U> &other);

  template <typename U>
  friend std::ostream &operator<<(std::ostream &l, const Matrix<U> &m);

  template <typename U> friend void operator>>(std::istream &l, Matrix<U> &m);

  template <typename U>
  friend matrix::Matrix<U> operator*(matrix::Matrix<U> &This, U other);

  template <typename U>
  friend matrix::Matrix<U> operator*(matrix::Matrix<U> &This, matrix::Matrix<U> &other);

  template <typename U>
  friend Matrix<U> operator+(Matrix<U> &This, Matrix<U> &other);

  template <typename U>
  friend Matrix<U> operator-(Matrix<U> &This, Matrix<U> &other);

private:
  std::vector<T> data;
  std::pair<int, int> dimensions;

public:
  static const Matrix<T> MZero;
  static const Matrix<T> MOne;
};

} // namespace matrix
