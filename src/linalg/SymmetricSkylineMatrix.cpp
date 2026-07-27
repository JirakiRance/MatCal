#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"

namespace MatCal::Linalg {

SymmetricSkylineMatrix::SymmetricSkylineMatrix(size_type size) {
    std::vector<size_type> first_columns(size);
    for (size_type i = 0; i < size; ++i) {
        first_columns[i] = i;
    }
    initialize_from_profile(std::move(first_columns));
}

SymmetricSkylineMatrix::SymmetricSkylineMatrix(std::vector<size_type> first_columns) {
    initialize_from_profile(std::move(first_columns));
}

SymmetricSkylineMatrix SymmetricSkylineMatrix::from_profile(std::vector<size_type> first_columns) {
    return SymmetricSkylineMatrix(std::move(first_columns));
}

SymmetricSkylineMatrix SymmetricSkylineMatrix::from_symmetric_positions(size_type size,
                                                                        const std::vector<index_pair>& positions) {
    std::vector<size_type> first_columns(size);
    for (size_type i = 0; i < size; ++i) {
        first_columns[i] = i;
    }

    for (const auto& [row, column] : positions) {
        if (row >= size || column >= size) {
            throw std::out_of_range("Skyline position out of range");
        }
        size_type high = std::max(row, column);
        size_type low = std::min(row, column);
        first_columns[high] = std::min(first_columns[high], low);
    }

    return SymmetricSkylineMatrix(std::move(first_columns));
}

bool SymmetricSkylineMatrix::stores(size_type row, size_type column) const {
    check_index(row);
    check_index(column);
    size_type high = std::max(row, column);
    size_type low = std::min(row, column);
    return low >= first_columns_[high];
}

double SymmetricSkylineMatrix::get(size_type row, size_type column) const {
    check_index(row);
    check_index(column);
    size_type high = std::max(row, column);
    size_type low = std::min(row, column);
    if (low < first_columns_[high]) {
        return 0.0;
    }
    return values_[storage_index(high, low)];
}

void SymmetricSkylineMatrix::set(size_type row, size_type column, double value) {
    check_index(row);
    check_index(column);
    size_type high = std::max(row, column);
    size_type low = std::min(row, column);
    if (low < first_columns_[high]) {
        throw std::out_of_range("Skyline set outside profile");
    }
    values_[storage_index(high, low)] = value;
}

void SymmetricSkylineMatrix::add(size_type row, size_type column, double value) {
    check_index(row);
    check_index(column);
    size_type high = std::max(row, column);
    size_type low = std::min(row, column);
    if (low < first_columns_[high]) {
        throw std::out_of_range("Skyline add outside profile");
    }
    values_[storage_index(high, low)] += value;
}

Vector SymmetricSkylineMatrix::multiply(const Vector& x) const {
    if (x.size() != size_) {
        throw std::invalid_argument("Skyline matvec dimension mismatch");
    }

    Vector result(size_);
    for (size_type row = 0; row < size_; ++row) {
        for (size_type column = first_columns_[row]; column <= row; ++column) {
            double value = values_[storage_index(row, column)];
            result[row] += value * x[column];
            if (column != row) {
                result[column] += value * x[row];
            }
        }
    }
    return result;
}

double SymmetricSkylineMatrix::matrix_scale() const noexcept {
    double scale = 0.0;
    for (double value : values_) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

bool SymmetricSkylineMatrix::all_finite() const noexcept {
    return std::all_of(values_.begin(), values_.end(), [](double value) {
        return std::isfinite(value);
    });
}

void SymmetricSkylineMatrix::initialize_from_profile(std::vector<size_type> first_columns) {
    size_ = first_columns.size();
    first_columns_ = std::move(first_columns);
    row_offsets_.assign(size_ + 1, 0);

    size_type storage = 0;
    for (size_type row = 0; row < size_; ++row) {
        if (first_columns_[row] > row) {
            throw std::invalid_argument("Skyline profile first column must be <= row");
        }
        row_offsets_[row] = storage;
        size_type width = checked_profile_width(row, first_columns_[row]);
        if (storage > std::numeric_limits<size_type>::max() - width) {
            throw std::length_error("Skyline storage size overflow");
        }
        storage += width;
    }
    row_offsets_[size_] = storage;
    values_.assign(storage, 0.0);
}

void SymmetricSkylineMatrix::check_index(size_type index) const {
    if (index >= size_) {
        throw std::out_of_range("Skyline index out of range");
    }
}

SymmetricSkylineMatrix::size_type SymmetricSkylineMatrix::storage_index(size_type row, size_type column) const {
    return row_offsets_[row] + (column - first_columns_[row]);
}

SymmetricSkylineMatrix::size_type SymmetricSkylineMatrix::checked_profile_width(size_type row,
                                                                                size_type first_column) const {
    return row - first_column + 1;
}

} // namespace MatCal::Linalg
