#ifndef LIBMATRIX_TEST_HPP
#define LIBMATRIX_TEST_HPP

#include "matrix_impl.hpp"

namespace test {

class MatrixTestCase {
public:
  MatrixTestCase() = delete;
  ~MatrixTestCase() = delete;

  static void run_tests();

  static int test_multadd_rows();

  static int test_normalize_row();

  static int test_swap_row();

  static int test_gaus_uniform();

  static int test_gaus_int();

  static int test_multiplication();
  
  static int test_multiplication2();

};

} // namespace test

#endif // LIBMATRIX_TEST_HPP
