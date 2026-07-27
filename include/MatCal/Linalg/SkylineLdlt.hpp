#ifndef MATCAL_LINALG_SKYLINE_LDLT_HPP
#define MATCAL_LINALG_SKYLINE_LDLT_HPP

#include <vector>

#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"

namespace MatCal::Linalg {

struct SkylineLdltFactorizationResult;

class SkylineLdltFactorization {
public:
    SkylineLdltFactorization() = default;

    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    std::size_t storage_size() const noexcept { return l_values_.size(); }
    double matrix_scale() const noexcept { return matrix_scale_; }
    double pivot_tolerance_used() const noexcept { return pivot_tolerance_used_; }
    double minimum_abs_pivot() const noexcept { return minimum_abs_pivot_; }

    double diagonal(std::size_t index) const;
    double lower(std::size_t row, std::size_t column) const;
    Vector multiply(const Vector& x) const;
    SolverResult solve(const Vector& b, const SolverOptions& options = {}) const;

private:
    friend struct SkylineLdltFactorizationResult;
    friend SkylineLdltFactorizationResult factorize_skyline_ldlt(const SymmetricSkylineMatrix&,
                                                                 const SolverOptions&);

    bool stores(std::size_t row, std::size_t column) const;
    std::size_t storage_index(std::size_t row, std::size_t column) const;

    std::size_t size_ = 0;
    std::vector<std::size_t> first_columns_;
    std::vector<std::size_t> row_offsets_;
    std::vector<double> l_values_;
    std::vector<double> diagonal_;
    double matrix_scale_ = 0.0;
    double pivot_tolerance_used_ = 0.0;
    double minimum_abs_pivot_ = 0.0;
    std::size_t factorization_operation_count_ = 0;
};

struct SkylineLdltFactorizationResult {
    SolverStatus status = SolverStatus::invalid_input;
    SkylineLdltFactorization factorization;
    SolverMetrics metrics;
    std::vector<SolverDiagnostic> diagnostics;
    std::string method = "skyline_ldlt";
    std::string implementation = "MatCal::Linalg skyline reference";

    bool success() const noexcept { return status == SolverStatus::success; }
};

SkylineLdltFactorizationResult factorize_skyline_ldlt(const SymmetricSkylineMatrix& matrix,
                                                       const SolverOptions& options = {});

SolverResult solve_skyline_ldlt(const SymmetricSkylineMatrix& matrix,
                                const Vector& b,
                                const SolverOptions& options = {});

} // namespace MatCal::Linalg

#endif
