#ifndef MATCAL_INTERPOLATION_LINEAR_INTERPOLATOR_HPP
#define MATCAL_INTERPOLATION_LINEAR_INTERPOLATOR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MatCal::Interpolation {

enum class ExtrapolationPolicy {
    reject,
    clamp,
    extrapolate
};

class LinearInterpolator {
public:
    LinearInterpolator() = default;

    LinearInterpolator(std::vector<double> xs,
                       std::vector<double> ys,
                       ExtrapolationPolicy policy = ExtrapolationPolicy::reject)
        : xs_(std::move(xs)),
          ys_(std::move(ys)),
          policy_(policy) {
        validate_nodes();
    }

    LinearInterpolator(const std::vector<std::pair<double, double>>& points,
                       ExtrapolationPolicy policy = ExtrapolationPolicy::reject)
        : policy_(policy) {
        xs_.reserve(points.size());
        ys_.reserve(points.size());
        for (const auto& [x, y] : points) {
            xs_.push_back(x);
            ys_.push_back(y);
        }
        validate_nodes();
    }

    LinearInterpolator(std::initializer_list<std::pair<double, double>> points,
                       ExtrapolationPolicy policy = ExtrapolationPolicy::reject)
        : LinearInterpolator(std::vector<std::pair<double, double>>(points), policy) {}

    double evaluate(double x) const {
        if (xs_.empty()) {
            throw std::logic_error("linear interpolator has no nodes");
        }
        check_finite(x, "linear interpolation point must be finite");

        if (x < xs_.front()) {
            return evaluate_outside_left(x);
        }
        if (x > xs_.back()) {
            return evaluate_outside_right(x);
        }
        if (x == xs_.back()) {
            return ys_.back();
        }

        const std::size_t i = interval_index(x);
        return interpolate(i, x);
    }

    const std::vector<double>& xs() const noexcept {
        return xs_;
    }

    const std::vector<double>& ys() const noexcept {
        return ys_;
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

    void validate_nodes() {
        if (xs_.size() != ys_.size()) {
            throw std::invalid_argument("linear interpolation xs and ys sizes must match");
        }
        if (xs_.size() < 2) {
            throw std::invalid_argument("linear interpolation requires at least two nodes");
        }
        for (std::size_t i = 0; i < xs_.size(); ++i) {
            check_finite(xs_[i], "linear interpolation x node must be finite");
            check_finite(ys_[i], "linear interpolation y node must be finite");
            if (i > 0 && !(xs_[i - 1] < xs_[i])) {
                throw std::invalid_argument("linear interpolation x nodes must be strictly increasing");
            }
        }
    }

    std::size_t interval_index(double x) const {
        auto upper = std::upper_bound(xs_.begin(), xs_.end(), x);
        const auto index = static_cast<std::size_t>(std::distance(xs_.begin(), upper));
        return index == 0 ? 0 : index - 1;
    }

    double interpolate(std::size_t i, double x) const {
        const double x0 = xs_[i];
        const double x1 = xs_[i + 1];
        const double y0 = ys_[i];
        const double y1 = ys_[i + 1];
        const double t = (x - x0) / (x1 - x0);
        const double value = y0 + t * (y1 - y0);
        check_finite(value, "linear interpolation produced a non-finite value");
        return value;
    }

    double evaluate_outside_left(double x) const {
        if (policy_ == ExtrapolationPolicy::reject) {
            throw std::out_of_range("linear interpolation point is left of the node range");
        }
        if (policy_ == ExtrapolationPolicy::clamp) {
            return ys_.front();
        }
        return interpolate(0, x);
    }

    double evaluate_outside_right(double x) const {
        if (policy_ == ExtrapolationPolicy::reject) {
            throw std::out_of_range("linear interpolation point is right of the node range");
        }
        if (policy_ == ExtrapolationPolicy::clamp) {
            return ys_.back();
        }
        return interpolate(xs_.size() - 2, x);
    }

    std::vector<double> xs_;
    std::vector<double> ys_;
    ExtrapolationPolicy policy_ = ExtrapolationPolicy::reject;
};

} // namespace MatCal::Interpolation

#endif
