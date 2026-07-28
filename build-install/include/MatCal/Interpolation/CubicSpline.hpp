#ifndef MATCAL_INTERPOLATION_CUBIC_SPLINE_HPP
#define MATCAL_INTERPOLATION_CUBIC_SPLINE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MatCal/Interpolation/LinearInterpolator.hpp"

namespace MatCal::Interpolation {

class CubicSpline {
public:
    CubicSpline() = default;

    CubicSpline(std::vector<double> xs,
                std::vector<double> ys,
                ExtrapolationPolicy policy = ExtrapolationPolicy::reject)
        : xs_(std::move(xs)),
          ys_(std::move(ys)),
          policy_(policy) {
        construct();
    }

    double evaluate(double x) const {
        if (xs_.empty()) {
            throw std::logic_error("cubic spline has no nodes");
        }
        check_finite(x, "cubic spline evaluation point must be finite");

        if (x < xs_.front()) {
            if (policy_ == ExtrapolationPolicy::reject) {
                throw std::out_of_range("cubic spline point is left of the node range");
            }
            if (policy_ == ExtrapolationPolicy::clamp) {
                return ys_.front();
            }
        }
        if (x > xs_.back()) {
            if (policy_ == ExtrapolationPolicy::reject) {
                throw std::out_of_range("cubic spline point is right of the node range");
            }
            if (policy_ == ExtrapolationPolicy::clamp) {
                return ys_.back();
            }
        }
        if (x == xs_.back()) {
            return ys_.back();
        }

        const std::size_t i = interval_index(x);
        const double h = xs_[i + 1] - xs_[i];
        const double a = (xs_[i + 1] - x) / h;
        const double b = (x - xs_[i]) / h;
        const double value = a * ys_[i] + b * ys_[i + 1] +
                             ((a * a * a - a) * h * h / 6.0) * second_derivatives_[i] +
                             ((b * b * b - b) * h * h / 6.0) * second_derivatives_[i + 1];
        check_finite(value, "cubic spline evaluation produced a non-finite value");
        return value;
    }

    double derivative(double x) const {
        if (xs_.empty()) {
            throw std::logic_error("cubic spline has no nodes");
        }
        check_finite(x, "cubic spline derivative point must be finite");
        if (x < xs_.front() || x > xs_.back()) {
            if (policy_ == ExtrapolationPolicy::reject) {
                throw std::out_of_range("cubic spline derivative point is outside the node range");
            }
            if (policy_ == ExtrapolationPolicy::clamp) {
                x = std::clamp(x, xs_.front(), xs_.back());
            }
        }
        if (x == xs_.back()) {
            x = xs_[xs_.size() - 2];
        }
        const std::size_t i = interval_index(x);
        const double h = xs_[i + 1] - xs_[i];
        const double a = (xs_[i + 1] - x) / h;
        const double b = (x - xs_[i]) / h;
        const double slope = (ys_[i + 1] - ys_[i]) / h;
        const double value = slope +
                             (h / 6.0) * ((1.0 - 3.0 * a * a) * second_derivatives_[i] +
                                          (3.0 * b * b - 1.0) * second_derivatives_[i + 1]);
        check_finite(value, "cubic spline derivative produced a non-finite value");
        return value;
    }

    const std::vector<double>& xs() const noexcept {
        return xs_;
    }

    const std::vector<double>& ys() const noexcept {
        return ys_;
    }

    const std::vector<double>& second_derivatives() const noexcept {
        return second_derivatives_;
    }

    ExtrapolationPolicy extrapolation_policy() const noexcept {
        return policy_;
    }

private:
    static void check_finite(double value, const char* message) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(message);
        }
    }

    void construct() {
        validate_nodes();
        second_derivatives_.assign(xs_.size(), 0.0);
        if (xs_.size() == 2) {
            return;
        }

        const std::size_t n = xs_.size();
        std::vector<double> lower(n, 0.0);
        std::vector<double> diag(n, 0.0);
        std::vector<double> upper(n, 0.0);
        std::vector<double> rhs(n, 0.0);

        diag[0] = 1.0;
        diag[n - 1] = 1.0;
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double h_left = xs_[i] - xs_[i - 1];
            const double h_right = xs_[i + 1] - xs_[i];
            lower[i] = h_left / 6.0;
            diag[i] = (h_left + h_right) / 3.0;
            upper[i] = h_right / 6.0;
            rhs[i] = (ys_[i + 1] - ys_[i]) / h_right -
                     (ys_[i] - ys_[i - 1]) / h_left;
            check_finite(rhs[i], "cubic spline rhs produced a non-finite value");
        }

        solve_tridiagonal(lower, diag, upper, rhs);
        second_derivatives_ = std::move(rhs);
        for (double value : second_derivatives_) {
            check_finite(value, "cubic spline second derivative is not finite");
        }
    }

    void validate_nodes() const {
        if (xs_.size() != ys_.size()) {
            throw std::invalid_argument("cubic spline xs and ys sizes must match");
        }
        if (xs_.size() < 2) {
            throw std::invalid_argument("cubic spline requires at least two nodes");
        }
        for (std::size_t i = 0; i < xs_.size(); ++i) {
            check_finite(xs_[i], "cubic spline x node must be finite");
            check_finite(ys_[i], "cubic spline y node must be finite");
            if (i > 0 && !(xs_[i - 1] < xs_[i])) {
                throw std::invalid_argument("cubic spline x nodes must be strictly increasing");
            }
        }
    }

    static void solve_tridiagonal(std::vector<double>& lower,
                                  std::vector<double>& diag,
                                  std::vector<double>& upper,
                                  std::vector<double>& rhs) {
        const std::size_t n = diag.size();
        for (std::size_t i = 1; i < n; ++i) {
            if (std::abs(diag[i - 1]) <= 1.0e-14) {
                throw std::runtime_error("cubic spline tridiagonal pivot is near zero");
            }
            const double factor = lower[i] / diag[i - 1];
            diag[i] -= factor * upper[i - 1];
            rhs[i] -= factor * rhs[i - 1];
            check_finite(diag[i], "cubic spline tridiagonal diagonal became non-finite");
            check_finite(rhs[i], "cubic spline tridiagonal rhs became non-finite");
        }

        if (std::abs(diag[n - 1]) <= 1.0e-14) {
            throw std::runtime_error("cubic spline tridiagonal final pivot is near zero");
        }
        rhs[n - 1] /= diag[n - 1];
        check_finite(rhs[n - 1], "cubic spline back substitution became non-finite");
        for (std::size_t offset = 1; offset < n; ++offset) {
            const std::size_t i = n - 1 - offset;
            if (std::abs(diag[i]) <= 1.0e-14) {
                throw std::runtime_error("cubic spline tridiagonal pivot is near zero");
            }
            rhs[i] = (rhs[i] - upper[i] * rhs[i + 1]) / diag[i];
            check_finite(rhs[i], "cubic spline back substitution became non-finite");
        }
    }

    std::size_t interval_index(double x) const {
        if (x <= xs_.front()) return 0;
        if (x >= xs_.back()) return xs_.size() - 2;
        auto upper = std::upper_bound(xs_.begin(), xs_.end(), x);
        const auto index = static_cast<std::size_t>(std::distance(xs_.begin(), upper));
        return index == 0 ? 0 : index - 1;
    }

    std::vector<double> xs_;
    std::vector<double> ys_;
    std::vector<double> second_derivatives_;
    ExtrapolationPolicy policy_ = ExtrapolationPolicy::reject;
};

} // namespace MatCal::Interpolation

#endif
