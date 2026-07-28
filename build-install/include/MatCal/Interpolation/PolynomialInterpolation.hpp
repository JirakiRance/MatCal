#ifndef MATCAL_INTERPOLATION_POLYNOMIAL_INTERPOLATION_HPP
#define MATCAL_INTERPOLATION_POLYNOMIAL_INTERPOLATION_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MatCal/Polynomial/Polynomial.hpp"

namespace MatCal::Interpolation {

struct DividedDifferenceResult {
    MatCal::Polynomial::Polynomial polynomial;
    std::vector<double> xs;
    std::vector<std::vector<double>> table;
};

namespace detail {

inline void validate_points(const std::vector<std::pair<double, double>>& points) {
    if (points.size() < 2) {
        throw std::invalid_argument("polynomial interpolation requires at least two points");
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (!std::isfinite(points[i].first) || !std::isfinite(points[i].second)) {
            throw std::invalid_argument("polynomial interpolation points must be finite");
        }
        for (std::size_t j = i + 1; j < points.size(); ++j) {
            if (points[i].first == points[j].first) {
                throw std::invalid_argument("polynomial interpolation duplicate x node");
            }
        }
    }
}

inline MatCal::Polynomial::Polynomial linear_factor(double root) {
    return MatCal::Polynomial::Polynomial::from_terms({{1, 1.0}, {0, -root}});
}

inline MatCal::Polynomial::Polynomial newton_basis(const std::vector<double>& xs, std::size_t count) {
    MatCal::Polynomial::Polynomial result({1.0});
    for (std::size_t i = 0; i < count; ++i) {
        result = result * linear_factor(xs[i]);
    }
    return result;
}

inline double factorial_scaled_denominator(double h, std::size_t n) {
    double result = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        result *= h * static_cast<double>(k);
        if (!std::isfinite(result)) {
            throw std::overflow_error("finite-difference interpolation denominator overflow");
        }
    }
    return result;
}

} // namespace detail

inline MatCal::Polynomial::Polynomial interpolate_lagrange(const std::vector<std::pair<double, double>>& points) {
    detail::validate_points(points);
    MatCal::Polynomial::Polynomial result;
    for (std::size_t i = 0; i < points.size(); ++i) {
        double coefficient = points[i].second;
        MatCal::Polynomial::Polynomial basis({1.0});
        for (std::size_t j = 0; j < points.size(); ++j) {
            if (i == j) {
                continue;
            }
            const double denominator = points[i].first - points[j].first;
            if (denominator == 0.0) {
                throw std::invalid_argument("lagrange denominator is zero");
            }
            coefficient /= denominator;
            basis = basis * detail::linear_factor(points[j].first);
        }
        result = result + basis * coefficient;
    }
    return result;
}

inline DividedDifferenceResult interpolate_newton_divided(const std::vector<std::pair<double, double>>& points) {
    detail::validate_points(points);
    const std::size_t n = points.size();
    DividedDifferenceResult result;
    result.xs.resize(n);
    result.table.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        result.xs[i] = points[i].first;
        result.table[i][0] = points[i].second;
    }
    for (std::size_t col = 1; col < n; ++col) {
        for (std::size_t row = col; row < n; ++row) {
            const double denominator = result.xs[row] - result.xs[row - col];
            if (denominator == 0.0) {
                throw std::invalid_argument("divided-difference denominator is zero");
            }
            result.table[row][col] = (result.table[row][col - 1] - result.table[row - 1][col - 1]) / denominator;
            if (!std::isfinite(result.table[row][col])) {
                throw std::runtime_error("divided-difference interpolation produced a non-finite coefficient");
            }
        }
    }
    result.polynomial = MatCal::Polynomial::Polynomial::from_terms({{0, result.table[0][0]}});
    for (std::size_t i = 1; i < n; ++i) {
        result.polynomial = result.polynomial + detail::newton_basis(result.xs, i) * result.table[i][i];
    }
    return result;
}

inline DividedDifferenceResult interpolate_newton_finite(double h, double x0, const std::vector<double>& y_values) {
    if (!std::isfinite(h) || h <= 0.0 || !std::isfinite(x0)) {
        throw std::invalid_argument("finite-difference interpolation requires finite x0 and positive h");
    }
    if (y_values.size() < 2) {
        throw std::invalid_argument("finite-difference interpolation requires at least two values");
    }
    for (double y : y_values) {
        if (!std::isfinite(y)) {
            throw std::invalid_argument("finite-difference interpolation y values must be finite");
        }
    }

    const std::size_t n = y_values.size();
    DividedDifferenceResult result;
    result.xs.resize(n);
    result.table.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        result.xs[i] = x0 + static_cast<double>(i) * h;
        result.table[i][0] = y_values[i];
    }
    for (std::size_t col = 1; col < n; ++col) {
        for (std::size_t row = col; row < n; ++row) {
            result.table[row][col] = result.table[row][col - 1] - result.table[row - 1][col - 1];
        }
    }

    result.polynomial = MatCal::Polynomial::Polynomial::from_terms({{0, result.table[0][0]}});
    for (std::size_t i = 1; i < n; ++i) {
        const double coefficient = result.table[i][i] / detail::factorial_scaled_denominator(h, i);
        result.polynomial = result.polynomial + detail::newton_basis(result.xs, i) * coefficient;
    }
    return result;
}

inline MatCal::Polynomial::Polynomial interpolate_hermite(const std::vector<double>& xs,
                                                          const std::vector<double>& ys,
                                                          const std::vector<double>& derivatives) {
    if (xs.size() != ys.size() || xs.size() != derivatives.size() || xs.size() < 2) {
        throw std::invalid_argument("Hermite interpolation sizes must match and contain at least two nodes");
    }
    std::vector<std::pair<double, double>> points;
    points.reserve(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        points.emplace_back(xs[i], ys[i]);
    }
    detail::validate_points(points);
    for (double derivative : derivatives) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument("Hermite derivative values must be finite");
        }
    }

    std::vector<MatCal::Polynomial::Polynomial> lagrange_basis(xs.size());
    std::vector<double> lagrange_derivatives(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        double denominator = 1.0;
        MatCal::Polynomial::Polynomial basis({1.0});
        for (std::size_t j = 0; j < xs.size(); ++j) {
            if (i == j) {
                continue;
            }
            denominator *= xs[i] - xs[j];
            basis = basis * detail::linear_factor(xs[j]);
        }
        lagrange_basis[i] = basis / denominator;
        lagrange_derivatives[i] = lagrange_basis[i].derivative().evaluate(xs[i]);
    }

    MatCal::Polynomial::Polynomial result;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        auto x_minus_xi = detail::linear_factor(xs[i]);
        auto li_sq = lagrange_basis[i] * lagrange_basis[i];
        auto hi = (MatCal::Polynomial::Polynomial({1.0}) - x_minus_xi * (2.0 * lagrange_derivatives[i])) * li_sq;
        auto htilde = x_minus_xi * li_sq;
        result = result + hi * ys[i] + htilde * derivatives[i];
    }
    return result;
}

} // namespace MatCal::Interpolation

#endif
