#ifndef MATCAL_LEASTSQUARES_LEASTSQUARES_HPP
#define MATCAL_LEASTSQUARES_LEASTSQUARES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/Vector.hpp"
#include "MatCal/Polynomial/Polynomial.hpp"

namespace MatCal::LeastSquares {

enum class LeastSquaresStatus {
    success,
    invalid_input,
    non_finite_input,
    rank_deficient,
    breakdown
};

enum class LeastSquaresReason {
    none,
    invalid_degree,
    size_mismatch,
    no_terms_selected,
    invalid_weight,
    non_finite_input,
    normal_equation_overflow,
    rank_deficient,
    linear_solve_failed
};

struct LeastSquaresDiagnostic {
    LeastSquaresStatus status = LeastSquaresStatus::success;
    LeastSquaresReason reason = LeastSquaresReason::none;
    std::string phase;
    std::size_t index = invalid_index();
    std::string message;

    static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }
};

struct LeastSquaresMetrics {
    std::size_t sample_count = 0;
    std::size_t term_count = 0;
    double residual_norm_inf = 0.0;
};

struct LeastSquaresResult {
    LeastSquaresStatus status = LeastSquaresStatus::invalid_input;
    std::vector<double> coefficients;
    std::vector<int> selected_degrees;
    MatCal::Polynomial::Polynomial polynomial;
    MatCal::Linalg::DenseMatrix normal_matrix;
    MatCal::Linalg::Vector rhs;
    LeastSquaresDiagnostic diagnostic;
    LeastSquaresMetrics metrics;

    bool success() const noexcept {
        return status == LeastSquaresStatus::success;
    }
};

namespace detail {

inline LeastSquaresDiagnostic diagnostic(LeastSquaresStatus status,
                                         LeastSquaresReason reason,
                                         std::string phase,
                                         std::string message,
                                         std::size_t index = LeastSquaresDiagnostic::invalid_index()) {
    return {status, reason, std::move(phase), index, std::move(message)};
}

inline LeastSquaresResult failure(LeastSquaresStatus status,
                                  LeastSquaresReason reason,
                                  std::string phase,
                                  std::string message,
                                  LeastSquaresMetrics metrics,
                                  std::size_t index = LeastSquaresDiagnostic::invalid_index()) {
    LeastSquaresResult result;
    result.status = status;
    result.metrics = metrics;
    result.diagnostic = diagnostic(status, reason, std::move(phase), std::move(message), index);
    return result;
}

inline bool finite_values(const std::vector<double>& values) noexcept {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

inline std::vector<int> all_degrees(int degree) {
    std::vector<int> degrees(static_cast<std::size_t>(degree) + 1);
    for (int i = 0; i <= degree; ++i) {
        degrees[static_cast<std::size_t>(i)] = i;
    }
    return degrees;
}

inline double pow_int(double x, int exponent) {
    double result = 1.0;
    for (int i = 0; i < exponent; ++i) {
        result *= x;
        if (!std::isfinite(result)) {
            throw std::overflow_error("least-squares power overflow");
        }
    }
    return result;
}

} // namespace detail

inline LeastSquaresResult fit_polynomial(const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         const std::vector<double>& weights,
                                         const std::vector<int>& selected_degrees) {
    LeastSquaresMetrics metrics;
    metrics.sample_count = x.size();
    metrics.term_count = selected_degrees.size();

    if (x.size() != y.size() || x.size() != weights.size()) {
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::size_mismatch,
                               "setup", "least-squares x, y, and weights sizes must match", metrics);
    }
    if (x.empty()) {
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::size_mismatch,
                               "setup", "least-squares requires at least one sample", metrics);
    }
    if (selected_degrees.empty()) {
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::no_terms_selected,
                               "setup", "least-squares requires at least one selected term", metrics);
    }
    if (!detail::finite_values(x) || !detail::finite_values(y) || !detail::finite_values(weights)) {
        return detail::failure(LeastSquaresStatus::non_finite_input, LeastSquaresReason::non_finite_input,
                               "setup", "least-squares inputs must be finite", metrics);
    }
    for (std::size_t i = 0; i < weights.size(); ++i) {
        if (weights[i] <= 0.0) {
            return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::invalid_weight,
                                   "setup", "least-squares weights must be positive", metrics, i);
        }
    }
    for (int degree : selected_degrees) {
        if (degree < 0) {
            return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::invalid_degree,
                                   "setup", "least-squares selected degrees must be non-negative", metrics);
        }
    }

    const std::size_t m = selected_degrees.size();
    MatCal::Linalg::DenseMatrix normal(m, m);
    MatCal::Linalg::Vector rhs(m);

    try {
        for (std::size_t sample = 0; sample < x.size(); ++sample) {
            std::vector<double> basis(m);
            for (std::size_t term = 0; term < m; ++term) {
                basis[term] = detail::pow_int(x[sample], selected_degrees[term]);
            }
            for (std::size_t row = 0; row < m; ++row) {
                rhs[row] += weights[sample] * y[sample] * basis[row];
                if (!std::isfinite(rhs[row])) {
                    throw std::overflow_error("least-squares rhs overflow");
                }
                for (std::size_t col = 0; col <= row; ++col) {
                    normal(row, col) += weights[sample] * basis[row] * basis[col];
                    if (!std::isfinite(normal(row, col))) {
                        throw std::overflow_error("least-squares normal matrix overflow");
                    }
                    normal(col, row) = normal(row, col);
                }
            }
        }
    } catch (const std::exception& e) {
        return detail::failure(LeastSquaresStatus::breakdown, LeastSquaresReason::normal_equation_overflow,
                               "normal_equations", e.what(), metrics);
    }

    auto linear_result = MatCal::Linalg::solve_dense_partial_pivot(normal, rhs);
    if (!linear_result.success()) {
        const bool rank = linear_result.status == MatCal::Linalg::SolverStatus::singular;
        auto status = rank ? LeastSquaresStatus::rank_deficient : LeastSquaresStatus::breakdown;
        auto reason = rank ? LeastSquaresReason::rank_deficient : LeastSquaresReason::linear_solve_failed;
        LeastSquaresResult result = detail::failure(status, reason, "linear_solve",
                                                    "least-squares normal equation solve failed", metrics);
        result.normal_matrix = normal;
        result.rhs = rhs;
        return result;
    }

    std::vector<double> coefficients(m);
    std::vector<MatCal::Polynomial::Polynomial::term> terms;
    terms.reserve(m);
    for (std::size_t i = 0; i < m; ++i) {
        coefficients[i] = linear_result.solution[i];
        terms.emplace_back(static_cast<std::size_t>(selected_degrees[i]), coefficients[i]);
    }

    LeastSquaresResult result;
    result.status = LeastSquaresStatus::success;
    result.coefficients = std::move(coefficients);
    result.selected_degrees = selected_degrees;
    result.polynomial = MatCal::Polynomial::Polynomial::from_terms(terms);
    result.normal_matrix = normal;
    result.rhs = rhs;
    result.metrics = metrics;
    result.diagnostic = detail::diagnostic(LeastSquaresStatus::success, LeastSquaresReason::none,
                                           "complete", "least-squares fit succeeded");

    for (std::size_t i = 0; i < x.size(); ++i) {
        result.metrics.residual_norm_inf =
            std::max(result.metrics.residual_norm_inf, std::abs(result.polynomial.evaluate(x[i]) - y[i]));
    }
    return result;
}

inline LeastSquaresResult fit_polynomial_degree(int degree,
                                                const std::vector<double>& x,
                                                const std::vector<double>& y,
                                                const std::vector<double>& weights) {
    if (degree < 0) {
        LeastSquaresMetrics metrics;
        metrics.sample_count = x.size();
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::invalid_degree,
                               "setup", "least-squares degree must be non-negative", metrics);
    }
    return fit_polynomial(x, y, weights, detail::all_degrees(degree));
}

inline LeastSquaresResult fit_polynomial_degree(int degree,
                                                const std::vector<double>& x,
                                                const std::vector<double>& y) {
    return fit_polynomial_degree(degree, x, y, std::vector<double>(x.size(), 1.0));
}

inline LeastSquaresResult fit_polynomial_selected(int degree,
                                                  const std::vector<double>& x,
                                                  const std::vector<double>& y,
                                                  const std::vector<double>& weights,
                                                  const std::vector<bool>& selects) {
    if (degree < 0) {
        LeastSquaresMetrics metrics;
        metrics.sample_count = x.size();
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::invalid_degree,
                               "setup", "least-squares degree must be non-negative", metrics);
    }
    if (selects.size() != static_cast<std::size_t>(degree) + 1) {
        LeastSquaresMetrics metrics;
        metrics.sample_count = x.size();
        return detail::failure(LeastSquaresStatus::invalid_input, LeastSquaresReason::size_mismatch,
                               "setup", "least-squares selects size must be degree + 1", metrics);
    }
    std::vector<int> selected;
    for (int i = 0; i <= degree; ++i) {
        if (selects[static_cast<std::size_t>(i)]) {
            selected.push_back(i);
        }
    }
    return fit_polynomial(x, y, weights, selected);
}

} // namespace MatCal::LeastSquares

#endif
