/************************************
  Class implementing the Thomas
  algorithm for solving tridiagonal
  matrix systems.
************************************/
#ifndef TRIDIAGONAL_MATRIX_SOLVER_H
#define TRIDIAGONAL_MATRIX_SOLVER_H

#include <complex>
#include <stdexcept>
#include <vector>

template <typename T> class TridiagonalMatrixSolver
{
public:
  TridiagonalMatrixSolver(const unsigned int n) : size(n)
  {
    //The factorization and solve loops assume at least two rows
    if (n < 2)
      throw std::invalid_argument("TridiagonalMatrixSolver needs at least 2 rows");

    //Allocate space
    upper_diag.resize(size);
    diag.resize(size);
    lower_diag.resize(size);
    upper_diag_prime.resize(size);
    rhs_prime.resize(size);
  }

  void initialize(const double *lower_diagonal, const double *diagonal,
                  const double *upper_diagonal)
  {
    unsigned int n = size;
    //Copy the data to internal vectors
    for (unsigned int i = 0; i < size; ++i)
    {
      upper_diag[i] = upper_diagonal[i];
      lower_diag[i] = lower_diagonal[i];
      diag[i] = diagonal[i];
    }

    //Do initial factorization
    upper_diag_prime[0] = upper_diag[0] / diag[0];
    for (unsigned int i = 1; i < n - 1; ++i)
    {
      upper_diag_prime[i] = upper_diag[i] / (diag[i] - lower_diag[i] * upper_diag_prime[i - 1]);
    }
  }

  void solve(const T *rhs, const unsigned int stride, T *x)
  {
    unsigned int n = size;
    rhs_prime[0] = rhs[0] / diag[0];

    for (unsigned int i = 1; i < n - 1; ++i)
      rhs_prime[i] = (rhs[stride * i] - lower_diag[i] * rhs_prime[i - 1]) /
                     (diag[i] - lower_diag[i] * upper_diag_prime[i - 1]);

    rhs_prime[n - 1] = (rhs[(n - 1) * stride] - lower_diag[(n - 1)] * rhs_prime[n - 2]) /
                       (diag[(n - 1)] - lower_diag[(n - 1)] * upper_diag_prime[n - 2]);

    x[(n - 1) * stride] = rhs_prime[n - 1];
    for (int i = n - 2; i >= 0; --i)
      x[i * stride] = rhs_prime[i] - upper_diag_prime[i] * x[(i + 1) * stride];
  }

private:
  unsigned int size;
  std::vector<double> upper_diag, diag, lower_diag;
  std::vector<double> upper_diag_prime;
  std::vector<T> rhs_prime;
};

#endif
