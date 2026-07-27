#ifndef MATCAL_LINALG_SYMMETRIC_SKYLINE_MATRIX_HPP
#define MATCAL_LINALG_SYMMETRIC_SKYLINE_MATRIX_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

class SymmetricSkylineMatrix {
public:
    using size_type = std::size_t;
    using index_pair = std::pair<size_type, size_type>;

    SymmetricSkylineMatrix() = default;
    explicit SymmetricSkylineMatrix(size_type size);
    explicit SymmetricSkylineMatrix(std::vector<size_type> first_columns);

    static SymmetricSkylineMatrix from_profile(std::vector<size_type> first_columns);
    static SymmetricSkylineMatrix from_symmetric_positions(size_type size,
                                                           const std::vector<index_pair>& positions);

    size_type size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    size_type storage_size() const noexcept { return values_.size(); }

    std::span<double> values() noexcept { return values_; }
    std::span<const double> values() const noexcept { return values_; }
    const std::vector<size_type>& first_columns() const noexcept { return first_columns_; }
    const std::vector<size_type>& row_offsets() const noexcept { return row_offsets_; }

    bool stores(size_type row, size_type column) const;
    double get(size_type row, size_type column) const;
    void set(size_type row, size_type column, double value);
    void add(size_type row, size_type column, double value);

    Vector multiply(const Vector& x) const;
    double matrix_scale() const noexcept;
    bool all_finite() const noexcept;

private:
    void initialize_from_profile(std::vector<size_type> first_columns);
    void check_index(size_type index) const;
    size_type storage_index(size_type row, size_type column) const;
    size_type checked_profile_width(size_type row, size_type first_column) const;

    size_type size_ = 0;
    std::vector<size_type> first_columns_;
    std::vector<size_type> row_offsets_;
    std::vector<double> values_;
};

} // namespace MatCal::Linalg

#endif
