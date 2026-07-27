#ifndef MATCAL_POLYNOMIAL_POLYNOMIAL_HPP
#define MATCAL_POLYNOMIAL_POLYNOMIAL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MatCal::Polynomial {

class Polynomial {
public:
    using size_type = std::size_t;
    using term = std::pair<size_type, double>;

    static constexpr double zero_tolerance = 1.0e-12;

    Polynomial() = default;

    explicit Polynomial(std::vector<double> coefficients)
        : coefficients_(std::move(coefficients)) {
        canonicalize();
    }

    Polynomial(std::initializer_list<double> coefficients)
        : coefficients_(coefficients) {
        canonicalize();
    }

    static Polynomial from_terms(const std::vector<term>& terms) {
        size_type degree = 0;
        for (const auto& [term_degree, coefficient] : terms) {
            check_finite(coefficient, "polynomial coefficient must be finite");
            degree = std::max(degree, term_degree);
        }

        if (terms.empty()) {
            return Polynomial();
        }
        if (degree == std::numeric_limits<size_type>::max()) {
            throw std::length_error("polynomial degree overflow");
        }

        std::vector<double> coefficients(degree + 1, 0.0);
        for (const auto& [term_degree, coefficient] : terms) {
            coefficients[term_degree] += coefficient;
            check_finite(coefficients[term_degree], "polynomial coefficient accumulation overflow");
        }
        return Polynomial(std::move(coefficients));
    }

    static Polynomial from_terms(std::initializer_list<term> terms) {
        return from_terms(std::vector<term>(terms));
    }

    bool is_zero() const noexcept {
        return coefficients_.empty();
    }

    size_type degree() const noexcept {
        return coefficients_.empty() ? 0 : coefficients_.size() - 1;
    }

    size_type term_count() const noexcept {
        return coefficients_.empty() ? 0 : coefficients_.size();
    }

    double coefficient(size_type degree) const noexcept {
        return degree < coefficients_.size() ? coefficients_[degree] : 0.0;
    }

    const std::vector<double>& coefficients() const noexcept {
        return coefficients_;
    }

    std::vector<term> terms_descending() const {
        std::vector<term> terms;
        for (size_type index = coefficients_.size(); index-- > 0;) {
            if (std::abs(coefficients_[index]) >= zero_tolerance) {
                terms.emplace_back(index, coefficients_[index]);
            }
        }
        return terms;
    }

    double evaluate(double x) const {
        check_finite(x, "polynomial evaluation point must be finite");
        double result = 0.0;
        for (size_type index = coefficients_.size(); index-- > 0;) {
            result = result * x + coefficients_[index];
            check_finite(result, "polynomial evaluation produced a non-finite value");
        }
        return result;
    }

    Polynomial derivative() const {
        if (coefficients_.size() <= 1) {
            return Polynomial();
        }

        std::vector<double> result(coefficients_.size() - 1, 0.0);
        for (size_type degree = 1; degree < coefficients_.size(); ++degree) {
            result[degree - 1] = coefficients_[degree] * static_cast<double>(degree);
            check_finite(result[degree - 1], "polynomial derivative coefficient overflow");
        }
        return Polynomial(std::move(result));
    }

    Polynomial integral(double constant = 0.0) const {
        check_finite(constant, "polynomial integral constant must be finite");
        if (coefficients_.size() == std::numeric_limits<size_type>::max()) {
            throw std::length_error("polynomial integral degree overflow");
        }

        std::vector<double> result(coefficients_.size() + 1, 0.0);
        result[0] = constant;
        for (size_type degree = 0; degree < coefficients_.size(); ++degree) {
            result[degree + 1] = coefficients_[degree] / static_cast<double>(degree + 1);
            check_finite(result[degree + 1], "polynomial integral coefficient overflow");
        }
        return Polynomial(std::move(result));
    }

    double definite_integral(double a, double b) const {
        check_finite(a, "polynomial integral lower bound must be finite");
        check_finite(b, "polynomial integral upper bound must be finite");
        Polynomial antiderivative = integral();
        return antiderivative.evaluate(b) - antiderivative.evaluate(a);
    }

    std::function<double(double)> to_function() const {
        Polynomial copy = *this;
        return [copy](double x) {
            return copy.evaluate(x);
        };
    }

    friend Polynomial operator+(const Polynomial& left, const Polynomial& right) {
        std::vector<double> result(std::max(left.coefficients_.size(), right.coefficients_.size()), 0.0);
        for (size_type i = 0; i < result.size(); ++i) {
            result[i] = left.coefficient(i) + right.coefficient(i);
            check_finite(result[i], "polynomial addition coefficient overflow");
        }
        return Polynomial(std::move(result));
    }

    friend Polynomial operator-(const Polynomial& left, const Polynomial& right) {
        std::vector<double> result(std::max(left.coefficients_.size(), right.coefficients_.size()), 0.0);
        for (size_type i = 0; i < result.size(); ++i) {
            result[i] = left.coefficient(i) - right.coefficient(i);
            check_finite(result[i], "polynomial subtraction coefficient overflow");
        }
        return Polynomial(std::move(result));
    }

    friend Polynomial operator*(const Polynomial& left, const Polynomial& right) {
        if (left.is_zero() || right.is_zero()) {
            return Polynomial();
        }
        if (left.degree() > std::numeric_limits<size_type>::max() - right.degree()) {
            throw std::length_error("polynomial multiplication degree overflow");
        }

        std::vector<double> result(left.degree() + right.degree() + 1, 0.0);
        for (size_type i = 0; i < left.coefficients_.size(); ++i) {
            for (size_type j = 0; j < right.coefficients_.size(); ++j) {
                result[i + j] += left.coefficients_[i] * right.coefficients_[j];
                check_finite(result[i + j], "polynomial multiplication coefficient overflow");
            }
        }
        return Polynomial(std::move(result));
    }

    friend Polynomial operator*(const Polynomial& polynomial, double scalar) {
        check_finite(scalar, "polynomial scalar must be finite");
        std::vector<double> result = polynomial.coefficients_;
        for (double& coefficient : result) {
            coefficient *= scalar;
            check_finite(coefficient, "polynomial scalar multiplication overflow");
        }
        return Polynomial(std::move(result));
    }

    friend Polynomial operator*(double scalar, const Polynomial& polynomial) {
        return polynomial * scalar;
    }

    friend Polynomial operator/(const Polynomial& polynomial, double scalar) {
        check_finite(scalar, "polynomial divisor must be finite");
        if (std::abs(scalar) < zero_tolerance) {
            throw std::invalid_argument("polynomial division by zero");
        }
        return polynomial * (1.0 / scalar);
    }

private:
    static void check_finite(double value, const char* message) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(message);
        }
    }

    void canonicalize() {
        for (double coefficient : coefficients_) {
            check_finite(coefficient, "polynomial coefficient must be finite");
        }
        while (!coefficients_.empty() && std::abs(coefficients_.back()) < zero_tolerance) {
            coefficients_.pop_back();
        }
    }

    std::vector<double> coefficients_;
};

} // namespace MatCal::Polynomial

#endif
