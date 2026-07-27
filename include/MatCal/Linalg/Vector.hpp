#ifndef MATCAL_LINALG_VECTOR_HPP
#define MATCAL_LINALG_VECTOR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MatCal::Linalg {

class Vector {
public:
    using value_type = double;
    using size_type = std::size_t;
    using iterator = std::vector<double>::iterator;
    using const_iterator = std::vector<double>::const_iterator;

    Vector() = default;
    explicit Vector(size_type size) : values_(size) {}
    Vector(size_type size, double value) : values_(size, value) {}
    Vector(std::initializer_list<double> values) : values_(values) {}
    explicit Vector(std::vector<double> values) : values_(std::move(values)) {}

    size_type size() const noexcept { return values_.size(); }
    bool empty() const noexcept { return values_.empty(); }

    double* data() noexcept { return values_.data(); }
    const double* data() const noexcept { return values_.data(); }

    std::span<double> span() noexcept { return values_; }
    std::span<const double> span() const noexcept { return values_; }

    double& operator[](size_type index) noexcept { return values_[index]; }
    const double& operator[](size_type index) const noexcept { return values_[index]; }

    double& at(size_type index) {
        if (index >= size()) {
            throw std::out_of_range("Vector index out of range");
        }
        return values_[index];
    }

    const double& at(size_type index) const {
        if (index >= size()) {
            throw std::out_of_range("Vector index out of range");
        }
        return values_[index];
    }

    iterator begin() noexcept { return values_.begin(); }
    iterator end() noexcept { return values_.end(); }
    const_iterator begin() const noexcept { return values_.begin(); }
    const_iterator end() const noexcept { return values_.end(); }
    const_iterator cbegin() const noexcept { return values_.cbegin(); }
    const_iterator cend() const noexcept { return values_.cend(); }

    void fill(double value) {
        std::fill(values_.begin(), values_.end(), value);
    }

    bool all_finite() const noexcept {
        return std::all_of(values_.begin(), values_.end(), [](double value) {
            return std::isfinite(value);
        });
    }

    double dot(const Vector& other) const {
        require_same_size(other, "Vector dot size mismatch");
        double result = 0.0;
        for (size_type i = 0; i < size(); ++i) {
            result += values_[i] * other.values_[i];
        }
        return result;
    }

    double norm1() const {
        double result = 0.0;
        for (double value : values_) {
            result += std::abs(value);
        }
        return result;
    }

    double norm2() const {
        double scale = 0.0;
        double sumsq = 1.0;
        bool saw_nan = false;

        for (double value : values_) {
            double abs_value = std::abs(value);
            if (std::isnan(abs_value)) {
                saw_nan = true;
            } else if (std::isinf(abs_value)) {
                return abs_value;
            } else if (abs_value != 0.0) {
                if (scale < abs_value) {
                    double ratio = scale / abs_value;
                    sumsq = 1.0 + sumsq * ratio * ratio;
                    scale = abs_value;
                } else {
                    double ratio = abs_value / scale;
                    sumsq += ratio * ratio;
                }
            }
        }

        if (saw_nan) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return scale == 0.0 ? 0.0 : scale * std::sqrt(sumsq);
    }

    double normInf() const {
        double result = 0.0;
        for (double value : values_) {
            result = std::max(result, std::abs(value));
        }
        return result;
    }

    void scale(double alpha) {
        for (double& value : values_) {
            value *= alpha;
        }
    }

    void axpy(double alpha, const Vector& x) {
        require_same_size(x, "Vector axpy size mismatch");
        for (size_type i = 0; i < size(); ++i) {
            values_[i] += alpha * x.values_[i];
        }
    }

    const std::vector<double>& values() const noexcept { return values_; }

private:
    void require_same_size(const Vector& other, const char* message) const {
        if (size() != other.size()) {
            throw std::invalid_argument(message);
        }
    }

    std::vector<double> values_;
};

inline bool all_finite(const Vector& vector) noexcept {
    return vector.all_finite();
}

} // namespace MatCal::Linalg

#endif
