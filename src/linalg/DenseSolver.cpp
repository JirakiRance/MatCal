#include "MatCal/Linalg/DenseSolver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace MatCal::Linalg {
namespace {

SolverResult make_result(SolverStatus status, std::string message) {
    SolverResult result;
    result.status = status;
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg reference";
    if (!message.empty()) {
        result.diagnostics.push_back({status, std::move(message)});
    }
    return result;
}

void swap_rows(std::vector<double>& values,
               DenseMatrix::size_type cols,
               DenseMatrix::size_type row_a,
               DenseMatrix::size_type row_b) {
    if (row_a == row_b) {
        return;
    }
    for (DenseMatrix::size_type c = 0; c < cols; ++c) {
        std::swap(values[row_a * cols + c], values[row_b * cols + c]);
    }
}

} // namespace

double residual_norm_inf(const DenseMatrix& a, const Vector& x, const Vector& b) {
    if (a.cols() != x.size() || a.rows() != b.size()) {
        throw std::invalid_argument("residual dimension mismatch");
    }

    double max_residual = 0.0;
    for (DenseMatrix::size_type r = 0; r < a.rows(); ++r) {
        double sum = 0.0;
        for (DenseMatrix::size_type c = 0; c < a.cols(); ++c) {
            sum += a(r, c) * x[c];
        }
        max_residual = std::max(max_residual, std::abs(sum - b[r]));
    }
    return max_residual;
}

SolverResult solve_dense_partial_pivot(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options) {
    if (!options.valid()) {
        return make_result(SolverStatus::invalid_input, "Invalid solver options");
    }
    if (a.rows() != a.cols()) {
        return make_result(SolverStatus::dimension_mismatch, "Dense solve requires a square matrix");
    }
    if (a.rows() != b.size()) {
        return make_result(SolverStatus::dimension_mismatch, "Right-hand side size does not match matrix rows");
    }
    if (!a.all_finite() || !b.all_finite()) {
        return make_result(SolverStatus::non_finite_input, "Dense solve inputs must be finite");
    }

    const auto n = a.rows();
    SolverResult result;
    result.status = SolverStatus::success;
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg reference";
    result.solution = Vector(n);

    if (n == 0) {
        result.metrics.residual_norm = 0.0;
        return result;
    }

    std::vector<double> work = a.values();
    std::vector<double> rhs = b.values();
    const double matrix_scale = std::max(a.normInf(), 1.0);
    const double pivot_tolerance = options.pivot_tolerance(matrix_scale);

    for (DenseMatrix::size_type k = 0; k < n; ++k) {
        DenseMatrix::size_type pivot_row = k;
        double pivot_abs = std::abs(work[k * n + k]);
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            double candidate = std::abs(work[r * n + k]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = r;
            }
        }

        if (pivot_abs <= pivot_tolerance) {
            return make_result(SolverStatus::singular, "Pivot is zero or below pivot tolerance");
        }

        swap_rows(work, n, k, pivot_row);
        std::swap(rhs[k], rhs[pivot_row]);

        double pivot = work[k * n + k];
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            double factor = work[r * n + k] / pivot;
            work[r * n + k] = 0.0;
            for (DenseMatrix::size_type c = k + 1; c < n; ++c) {
                work[r * n + c] -= factor * work[k * n + c];
            }
            rhs[r] -= factor * rhs[k];
        }
        result.metrics.iterations = k + 1;
    }

    for (DenseMatrix::size_type i = n; i-- > 0;) {
        double sum = rhs[i];
        for (DenseMatrix::size_type c = i + 1; c < n; ++c) {
            sum -= work[i * n + c] * result.solution[c];
        }
        double diag = work[i * n + i];
        if (std::abs(diag) <= pivot_tolerance) {
            return make_result(SolverStatus::singular, "Back substitution pivot is below tolerance");
        }
        result.solution[i] = sum / diag;
    }

    result.metrics.residual_norm = residual_norm_inf(a, result.solution, b);
    if (!std::isfinite(result.metrics.residual_norm)) {
        return make_result(SolverStatus::breakdown, "Residual is not finite");
    }

    return result;
}

} // namespace MatCal::Linalg
