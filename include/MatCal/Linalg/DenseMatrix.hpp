#ifndef MATCAL_LINALG_DENSE_MATRIX_HPP
#define MATCAL_LINALG_DENSE_MATRIX_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

class DenseMatrix {
public:
    using value_type = double;
    using size_type = std::size_t;

    DenseMatrix() = default;

    DenseMatrix(size_type rows, size_type cols)
        : rows_(rows), cols_(cols), values_(checked_element_count(rows, cols)) {}

    DenseMatrix(size_type rows, size_type cols, double value)
        : rows_(rows), cols_(cols), values_(checked_element_count(rows, cols), value) {}

    DenseMatrix(std::initializer_list<std::initializer_list<double>> rows) {
        rows_ = rows.size();
        cols_ = rows_ == 0 ? 0 : rows.begin()->size();
        values_.reserve(checked_element_count(rows_, cols_));

        for (const auto& row : rows) {
            if (row.size() != cols_) {
                throw std::invalid_argument("DenseMatrix initializer rows must have equal length");
            }
            values_.insert(values_.end(), row.begin(), row.end());
        }
    }

    static DenseMatrix from_row_major(size_type rows, size_type cols, std::vector<double> values) {
        if (values.size() != checked_element_count(rows, cols)) {
            throw std::invalid_argument("DenseMatrix row-major data size mismatch");
        }
        DenseMatrix result;
        result.rows_ = rows;
        result.cols_ = cols;
        result.values_ = std::move(values);
        return result;
    }

    static DenseMatrix identity(size_type n) {
        DenseMatrix result(n, n);
        for (size_type i = 0; i < n; ++i) {
            result(i, i) = 1.0;
        }
        return result;
    }

    size_type rows() const noexcept { return rows_; }
    size_type cols() const noexcept { return cols_; }
    size_type size() const noexcept { return values_.size(); }
    bool empty() const noexcept { return values_.empty(); }

    double* data() noexcept { return values_.data(); }
    const double* data() const noexcept { return values_.data(); }

    std::span<double> span() noexcept { return values_; }
    std::span<const double> span() const noexcept { return values_; }

    std::span<double> row(size_type row_index) {
        check_row(row_index);
        if (cols_ == 0) {
            return std::span<double>();
        }
        return std::span<double>(values_.data() + row_index * cols_, cols_);
    }

    std::span<const double> row(size_type row_index) const {
        check_row(row_index);
        if (cols_ == 0) {
            return std::span<const double>();
        }
        return std::span<const double>(values_.data() + row_index * cols_, cols_);
    }

    double& operator()(size_type row_index, size_type col_index) noexcept {
        return values_[row_index * cols_ + col_index];
    }

    const double& operator()(size_type row_index, size_type col_index) const noexcept {
        return values_[row_index * cols_ + col_index];
    }

    double& at(size_type row_index, size_type col_index) {
        check_indices(row_index, col_index);
        return (*this)(row_index, col_index);
    }

    const double& at(size_type row_index, size_type col_index) const {
        check_indices(row_index, col_index);
        return (*this)(row_index, col_index);
    }

    void fill(double value) {
        std::fill(values_.begin(), values_.end(), value);
    }

    bool all_finite() const noexcept {
        return std::all_of(values_.begin(), values_.end(), [](double value) {
            return std::isfinite(value);
        });
    }

    DenseMatrix transpose() const {
        DenseMatrix result(cols_, rows_);
        for (size_type r = 0; r < rows_; ++r) {
            for (size_type c = 0; c < cols_; ++c) {
                result(c, r) = (*this)(r, c);
            }
        }
        return result;
    }

    Vector multiply(const Vector& vector) const {
        if (cols_ != vector.size()) {
            throw std::invalid_argument("DenseMatrix matvec dimension mismatch");
        }
        Vector result(rows_);
        for (size_type r = 0; r < rows_; ++r) {
            double sum = 0.0;
            for (size_type c = 0; c < cols_; ++c) {
                sum += (*this)(r, c) * vector[c];
            }
            result[r] = sum;
        }
        return result;
    }

    DenseMatrix multiply(const DenseMatrix& other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("DenseMatrix matmul dimension mismatch");
        }
        DenseMatrix result(rows_, other.cols_);
        for (size_type r = 0; r < rows_; ++r) {
            for (size_type k = 0; k < cols_; ++k) {
                double left = (*this)(r, k);
                for (size_type c = 0; c < other.cols_; ++c) {
                    result(r, c) += left * other(k, c);
                }
            }
        }
        return result;
    }

    double normInf() const {
        double max_row_sum = 0.0;
        for (size_type r = 0; r < rows_; ++r) {
            double row_sum = 0.0;
            for (size_type c = 0; c < cols_; ++c) {
                row_sum += std::abs((*this)(r, c));
            }
            max_row_sum = std::max(max_row_sum, row_sum);
        }
        return max_row_sum;
    }

    const std::vector<double>& values() const noexcept { return values_; }

    static size_type checked_element_count(size_type rows, size_type cols) {
        if (rows != 0 && cols > std::numeric_limits<size_type>::max() / rows) {
            throw std::length_error("DenseMatrix size overflow");
        }
        return rows * cols;
    }

private:
    void check_row(size_type row_index) const {
        if (row_index >= rows_) {
            throw std::out_of_range("DenseMatrix row index out of range");
        }
    }

    void check_indices(size_type row_index, size_type col_index) const {
        if (row_index >= rows_ || col_index >= cols_) {
            throw std::out_of_range("DenseMatrix index out of range");
        }
    }

    size_type rows_ = 0;
    size_type cols_ = 0;
    std::vector<double> values_;
};

inline bool all_finite(const DenseMatrix& matrix) noexcept {
    return matrix.all_finite();
}

} // namespace MatCal::Linalg

#endif
