#include "test.hpp"

int test::MatrixTestCase::test_multadd_rows() {
  matrix::Matrix<double> input(
      {{1.32, 2.10, 3.02}, {4.14, 1.1, 2.20}, {0.10, 0.9, 0.75}});

  matrix::Matrix<double> exptected(
      {{1.32, 2.10, 3.02}, {8.38, 3.1, 5.15}, {01.0, 0.9, 0.75}});

  input.multadd_rows(1, 2, -2);
  return input == exptected;
}

int test::MatrixTestCase::test_normalize_row() {
  matrix::Matrix<double> input({{1, 2, 3}, {4, 1, 2}, {1, 1, 1}});

  matrix::Matrix<double> exptected({{1, 2, 3}, {-8, -2, -4}, {1, 1, 1}});

  input.normalize_row(1, -2);

  return input == exptected;
}

int test::MatrixTestCase::test_swap_row() {
  matrix::Matrix<double> input({{1, 2, 3}, {4, 1, 2}, {0, 0, 0}});

  matrix::Matrix<double> exptected({{1, 2, 3}, {0, 0, 0}, {4, 1, 2}});

  input.swap_row(1, 2);
  return input == exptected;
}

int test::MatrixTestCase::test_gaus_uniform() {
  matrix::Matrix<double> input({{1.2, 4.135, 2.613, 31.18027},
                                {-3.851, -2.719, 16.172, -51, 40619},
                                {2.314, 23.187, 5.617, -15.98744}});

  matrix::Matrix<double> reduced = input.reduce();

  matrix::Matrix<double> exptected({{1, 0, 0, 34.6290859876895},
                                    {0, 1, 0, -5.16269041395079},
                                    {0, 0, 1, 4.19942276175242}});

//  std::cout << input;
//  std::cout << "====\n";

  return reduced == exptected;
}

int test::MatrixTestCase::test_gaus_int() {
  matrix::Matrix<double> input(
      {{2.0, 2, 1.0, 14}, {1, 1, -1, -11}, {4, 2, 3, 44}});

  matrix::Matrix<double> reduced = input.reduce();

  matrix::Matrix<double> exptected(
      {{1, 0, 0, 3}, {0, 1, 0, -2}, {0, 0, 1, 12}});

//  std::cout << input;
//  std::cout << "====\n";

  return reduced == exptected;
}

int test::MatrixTestCase::test_multiplication() {
  matrix::Matrix<int> x({{1, 2}, {1, 1}, {1, 4}});

  matrix::Matrix<int> y({{0, 2, 0, 1}, {2, 3, -1, 1}});

  matrix::Matrix<int> res = x * y;

  matrix::Matrix<int> expected_res(
      {{4, 8, 2, -1}, {2, 5, 1, 0}, {8, 14, 4, -3}});

  return res == expected_res;
}

int test::MatrixTestCase::test_multiplication2() {
  matrix::Matrix<int> x({{3,  5}});

  matrix::Matrix<int> y({{-1,  3,  -5}, {2, 10, -1}});

  matrix::Matrix<int> res = x * y;

  matrix::Matrix<int> expected_res(
      {{7,  59,  -20}});

  return res == expected_res;
}

void test::MatrixTestCase::run_tests() {
  int err;
  err = test_swap_row();
  std::cout << "test_swap_row " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_multadd_rows();
  std::cout << "test_multadd_rows " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_normalize_row();
  std::cout << "test_normalize_row " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_gaus_int();
  std::cout << "test_gaus_int " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_gaus_uniform();
  std::cout << "test_gaus_uniform " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_multiplication();
  std::cout << "test_multiplication " << (err ? "PASSED" : "FAILED") << "\n";

  err = test_multiplication2();
  std::cout << "test_multiplication2 " << (err ? "PASSED" : "FAILED") << "\n";
}

template class matrix::Matrix<double>;
