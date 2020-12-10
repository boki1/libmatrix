#ifndef LIBMATRIX_MATRIX_IMPL_HPP
#define LIBMATRIX_MATRIX_IMPL_HPP

#include <matrix.hpp>

namespace matrix {

template <typename T> Matrix<T>::Matrix() { this->set_empty(); }

template <typename T> Matrix<T>::Matrix(int d, int x, int y) {
  this->set_empty();
  this->set_dimensions(x, y);
  this->set_diagonal(d);
}

template <typename T> Matrix<T>::Matrix(std::vector<std::vector<T>> data) {
  this->set_empty();
  this->set_dimensions(data[0].size(), data.size());
  for (int row = 0; row < data.size(); ++row) {
    for (int col = 0; col < data[0].size(); ++col) {
      this->store_at(col, row, data[row][col]);
    }
  }
}

template <typename T>
Matrix<T> &Matrix<T>::normalize_row(int row, double coeff) {
  T temp;
  for (int col = 0; col < this->x_dim(); ++col) {
    temp = this->at(col, row) / coeff;
    this->store_at(col, row, temp);
  }
  return *this;
}

template <typename T> T *Matrix<T>::get_row(int i) const {
  T *_row = new T[x_dim()];
  for (int j = 0; j < x_dim(); ++j) {
    _row[j] = this->at(j, i);
  }
  return _row;
}

template <typename T> Matrix<T> &Matrix<T>::swap_row(int i, int j) {
  T *first = get_row(i), *second = get_row(j);
  for (int h = 0; h < x_dim(); ++h) {
    this->store_at(h, i, second[h]);
    this->store_at(h, j, first[h]);
  }

  delete[] first;
  delete[] second;
  return *this;
}
template <typename T> Matrix<T> &Matrix<T>::save_row(int i, T *values) {
  for (int j = 0; j < x_dim(); ++j) {
    this->store_at(j, i, values[j]);
  }

  return *this;
}
// `i` = `j` * val + `i` (vector version)
template <typename T> Matrix<T> &Matrix<T>::multadd_rows(int i, int j, T val) {
  T *ii = get_row(i), *jj = get_row(j);

  T *actual = new T[x_dim()];
  for (int h = 0; h < x_dim(); ++h) {
    actual[h] = jj[h] * val + ii[h];
  }

  this->save_row(i, actual);

  delete[] ii;
  delete[] jj;

  return *this;
}

// Returns the row echelon form of a matrix
template <typename T> Matrix<T> Matrix<T>::reduce() const {
  Matrix<T> nmatrix = *this;

  T pivot = 0;
  T val;
  int i;
  for (int r = 0; r < nmatrix.y_dim(); ++r) {
    // for (int r = 0; r < 1; ++r) {
    if (x_dim() - 1 <= pivot) {
      goto end;
    }

    i = r;
    while (nmatrix.at(pivot, i) == 0) {
      i++;
      if (y_dim() == i) {
        i = r;
        pivot++;
        if (x_dim() - 1 == pivot) {
          goto end;
        }
      }
    }

    if (r != i) {
      nmatrix.swap_row(i, r);
    }
    double norm_coeff = nmatrix.at(pivot, r);
    // std::cout << "nomalizing ...\n" << nmatrix << "row " << r <<
    // " with " << norm_coeff << "\n";
    nmatrix.normalize_row(r, nmatrix.at(pivot, r));

    // std::cout << "after normalizing\n" << nmatrix;

    // std::cout << "r = " << r << '\n';
    for (i = 0; i < y_dim(); ++i) {
      if (i != r) {
        val = -(nmatrix.at(pivot, i));
        // std::cout << "multiplying row " << i << "
        // with " << val << '\n';
        nmatrix.multadd_rows(i, r, val);
      }
    }
    pivot++;
  }
end:
  return nmatrix;
}

template <typename U> Matrix<U> operator+(Matrix<U> &This, Matrix<U> &other) {
  if (!This.same_dimensions_as(other))
    return Matrix<U>::MZero;

  Matrix<U> sum;
  sum->set_dimensions(This.x_dim(), This.y_dim());
  U val;
  for (int y = 0; y < This.y_dim(); ++y) {
    for (int x = 0; x < This.x_dim(); ++x) {
      val = This.at(x, y) + other.at(x, y);
      sum->store_at(x, y, val);
    }
  }

  return sum;
}

template <typename U>
std::ostream &operator<<(std::ostream &l, const Matrix<U> &m) {
  for (int y = 0; y < m.dimensions.second; ++y) {
    for (int x = 0; x < m.dimensions.first; ++x) {
      l << m.at(x, y) << ' ';
    }
    l << '\n';
  }
  return l;
}

template <typename U> void operator>>(std::istream &l, Matrix<U> &m) {
  int _y, _x;
  l >> _y >> _x;

  m.set_dimensions(_x, _y);
  U curr_val;

  for (int i = 0; i < _y; ++i) {
    for (int j = 0; j < _x; ++j) {
      l >> curr_val;
      m.store_at(j, i, curr_val);
    }
  }
}

template <typename U>
int operator==(matrix::Matrix<U> &This, matrix::Matrix<U> &other) {
  if (This.same_dimensions_as(other)) {
    for (int y = 0; y < This.dimensions.first; ++y) {
      for (int x = 0; x < other.dimensions.second; ++x) {
        U This_value = This.index(x, y);
        U other_value = other.index(x, y);
        if (This_value != other_value)
          return 0;
      }
    }
    return 1;
  }
  return 0;
}

template <typename U>
matrix::Matrix<U> operator*(matrix::Matrix<U> &This, U other) {
  U val;
  Matrix<U> res;
  for (int y = 0; y < This.y_dim(); ++y) {
    for (int x = 0; x < This.x_dim(); ++x) {
      val = This.at(x, y) * other;
      res.store_at(x, y, val);
    }
  }

  return res;
}

template <typename U> Matrix<U> operator-(Matrix<U> &This, Matrix<U> &other) {
  if (!This.same_dimensions_as(other))
    return Matrix<U>::MZero;

  Matrix<U> diff;
  diff->set_dimensions(This.x_dim(), This.y_dim());

  other = other * -1;
  *diff = This + other;
  other = other * -1;

  return diff;
}

// if C = AB for an n × m matrix A and an m × p matrix B, then C is an n × p
// matrix with entries
template <typename U>
matrix::Matrix<U> operator*(matrix::Matrix<U> &A, matrix::Matrix<U> &B) {
  auto dim_A = A.get_dimensions(), dim_B = B.get_dimensions();
  auto _m1 = dim_A.first;
  auto _m2 = dim_B.second;
  if (_m1 != _m2)
    throw matrix::UndefinedBahaviourMatrixException(
        "cannot multiply given matrices");

  Matrix<U> C;
  C.set_dimensions(dim_B.first, dim_A.second);
  int m, n, p;
  m = A.x_dim();
  n = A.y_dim();
  p = B.x_dim();
  U sum;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      sum = 0;
      for (int k = 0; k < m; ++k) {
        sum += A.at(k, i) * B.at(j, k);
      }
      C.store_at(i, j, sum);
    }
  }
  return C;
}

template <typename T> const Matrix<T> Matrix<T>::MZero = Matrix<T>(0);
template <typename T> const Matrix<T> Matrix<T>::MOne = Matrix<T>(1);

} // namespace matrix

#endif // LIBMATRIX_MATRIX_IMPL_HPP
