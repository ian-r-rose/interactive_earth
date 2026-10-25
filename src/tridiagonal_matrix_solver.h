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

/* Many tridiagonal systems of the same size, solved together: in the
   spectral solves there is one system per Fourier mode, running down a
   column of the (row, mode) array. Sweeping row by row with the systems
   innermost reads memory contiguously and vectorizes, unlike solving one
   column at a time.

   The matrices are real and factored once, keeping the reciprocals of the
   pivots, so a solve only multiplies. Each system is solved the same way
   no matter how the systems are split between calls (or threads). */
class BatchedTridiagonalSolver
{
public:
  BatchedTridiagonalSolver(int n_rows, int n_systems) : n_rows(n_rows), n_systems(n_systems)
  {
    //The factorization and solve loops assume at least two rows
    if (n_rows < 2)
      throw std::invalid_argument("BatchedTridiagonalSolver needs at least 2 rows");
    if (n_systems < 1)
      throw std::invalid_argument("BatchedTridiagonalSolver needs at least 1 system");

    lower.resize(n_rows * n_systems);
    inv_pivot.resize(n_rows * n_systems);
    upper_prime.resize(n_rows * n_systems);
  }

  //Set and factor the matrix of system s, given its diagonals (lower[0] and
  //upper[n_rows - 1] are unused)
  void initialize(int s, const double *lower_diagonal, const double *diagonal,
                  const double *upper_diagonal)
  {
    double previous_upper_prime = 0.0;
    for (int i = 0; i < n_rows; ++i)
    {
      const int k = i * n_systems + s;
      lower[k] = (i > 0 ? lower_diagonal[i] : 0.0);
      inv_pivot[k] = 1.0 / (diagonal[i] - lower[k] * previous_upper_prime);
      upper_prime[k] = (i < n_rows - 1 ? upper_diagonal[i] * inv_pivot[k] : 0.0);
      previous_upper_prime = upper_prime[k];
    }
  }

  //Solve systems first to last - 1, scaling the solution by scale. Element
  //(row i, system s) of rhs and x is at i * stride + s. x may be rhs.
  void solve(const std::complex<double> *rhs, std::complex<double> *x, int stride, int first,
             int last, double scale = 1.0) const
  {
    //Forward elimination
    for (int s = first; s < last; ++s)
      x[s] = scale * rhs[s] * inv_pivot[s];
    for (int i = 1; i < n_rows; ++i)
    {
      const std::complex<double> *b = rhs + i * stride;
      std::complex<double> *y = x + i * stride;
      const std::complex<double> *y_prev = y - stride;
      const double *l = &lower[i * n_systems], *p = &inv_pivot[i * n_systems];
      for (int s = first; s < last; ++s)
        y[s] = (scale * b[s] - l[s] * y_prev[s]) * p[s];
    }

    //Back substitution
    for (int i = n_rows - 2; i >= 0; --i)
    {
      std::complex<double> *y = x + i * stride;
      const std::complex<double> *y_next = y + stride;
      const double *u = &upper_prime[i * n_systems];
      for (int s = first; s < last; ++s)
        y[s] -= u[s] * y_next[s];
    }
  }

private:
  int n_rows, n_systems;
  //Indexed by (row i, system s) at i * n_systems + s
  std::vector<double> lower, inv_pivot, upper_prime;
};

#endif
