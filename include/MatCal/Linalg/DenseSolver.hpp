#ifndef MATCAL_LINALG_DENSE_SOLVER_HPP
#define MATCAL_LINALG_DENSE_SOLVER_HPP

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/Vector.hpp"

#include <vector>

namespace MatCal::Linalg {

double residual_norm_inf(const DenseMatrix& a, const Vector& x, const Vector& b);

struct MatrixSolverResult {
    SolverStatus status = SolverStatus::invalid_input;
    DenseMatrix solution;
    SolverMetrics metrics;
    std::vector<SolverDiagnostic> diagnostics;
    std::string method;
    std::string implementation;

    bool success() const noexcept {
        return status == SolverStatus::success;
    }
};

class PivotedLuFactorization {
public:
    using size_type = DenseMatrix::size_type;

    PivotedLuFactorization() = default;

    size_type size() const noexcept;
    bool valid() const noexcept;
    const DenseMatrix& lu() const noexcept;
    const DenseMatrix& original_matrix() const noexcept;
    const std::vector<size_type>& row_permutation() const noexcept;
    int permutation_sign() const noexcept;
    const SolverMetrics& factorization_metrics() const noexcept;

    double determinant() const;
    SolverResult solve(const Vector& rhs, const SolverOptions& options = {}) const;
    MatrixSolverResult solve(const DenseMatrix& rhs, const SolverOptions& options = {}) const;

private:
    friend struct PivotedLuBuilder;

    DenseMatrix original_;
    DenseMatrix lu_;
    std::vector<size_type> permutation_;
    int permutation_sign_ = 1;
    SolverMetrics metrics_;
    bool valid_ = false;
};

struct PivotedLuFactorizationResult {
    SolverStatus status = SolverStatus::invalid_input;
    PivotedLuFactorization factorization;
    SolverMetrics metrics;
    std::vector<SolverDiagnostic> diagnostics;
    std::string method;
    std::string implementation;

    bool success() const noexcept {
        return status == SolverStatus::success;
    }
};

PivotedLuFactorizationResult factorize_dense_partial_pivot(const DenseMatrix& a,
                                                           const SolverOptions& options = {});

SolverResult solve_dense_partial_pivot(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options = {});

} // namespace MatCal::Linalg

#endif
