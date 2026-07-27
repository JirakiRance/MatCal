#include <limits>
#include <vector>

#include "MatCal/Linalg/Vector.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

int main() {
    using MatCal::Linalg::Vector;

    Vector empty;
    expect_true(empty.empty() && empty.size() == 0, "default vector is empty");
    expect_near(empty.norm1(), 0.0, 0.0, "empty norm1");
    expect_near(empty.norm2(), 0.0, 0.0, "empty norm2");
    expect_near(empty.normInf(), 0.0, 0.0, "empty normInf");

    Vector sized(3);
    sized.fill(2.0);
    expect_true(sized.size() == 3, "size constructor");
    expect_near(sized[1], 2.0, 1e-12, "operator[] access");
    sized.at(1) = 4.0;
    expect_near(sized.at(1), 4.0, 1e-12, "checked at access");
    expect_throw([&] { (void)sized.at(3); }, "Vector at out of range");

    Vector listed{1.0, -2.0, 2.0};
    expect_near(listed.dot(Vector{3.0, 4.0, 5.0}), 5.0, 1e-12, "dot");
    expect_near(listed.norm1(), 5.0, 1e-12, "norm1");
    expect_near(listed.norm2(), 3.0, 1e-12, "norm2");
    expect_near(listed.normInf(), 2.0, 1e-12, "normInf");

    Vector y{1.0, 1.0, 1.0};
    y.axpy(2.0, listed);
    expect_near(y[0], 3.0, 1e-12, "axpy x0");
    expect_near(y[1], -3.0, 1e-12, "axpy x1");
    y.scale(0.5);
    expect_near(y[0], 1.5, 1e-12, "scale");

    std::vector<double> raw{4.0, 5.0};
    Vector from_vector(raw);
    expect_near(from_vector[1], 5.0, 1e-12, "explicit std::vector construction");
    from_vector.span()[0] = 6.0;
    expect_near(from_vector[0], 6.0, 1e-12, "writable span");
    std::span<const double> const_span = static_cast<const Vector&>(from_vector).span();
    expect_true(const_span.size() == 2, "const span");

    double iter_sum = 0.0;
    for (double value : listed) {
        iter_sum += value;
    }
    expect_near(iter_sum, 1.0, 1e-12, "begin end iteration");

    expect_throw([&] { (void)listed.dot(Vector{1.0}); }, "dot size mismatch");
    expect_throw([&] { listed.axpy(1.0, Vector{1.0}); }, "axpy size mismatch");

    Vector huge{1e200, 1e200};
    expect_true(std::isfinite(huge.norm2()), "stable norm2 avoids overflow");
    expect_relative(huge.norm2(), std::sqrt(2.0) * 1e200, 1e-12, "stable norm2 value");

    Vector overflow_norm{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    expect_true(std::isinf(overflow_norm.norm1()), "norm1 overflow is reported as Inf");
    expect_true(std::isinf(overflow_norm.dot(overflow_norm)), "dot overflow is reported as Inf");

    Vector non_finite{1.0, std::numeric_limits<double>::infinity()};
    expect_true(!non_finite.all_finite(), "finite check detects Inf");
    expect_true(std::isinf(non_finite.normInf()), "normInf reports Inf");
    Vector nan_values{std::numeric_limits<double>::quiet_NaN()};
    expect_true(!nan_values.all_finite(), "finite check detects NaN");
    expect_true(std::isnan(nan_values.normInf()), "normInf reports NaN");

    Vector overflow_scale{std::numeric_limits<double>::max()};
    overflow_scale.scale(2.0);
    expect_true(!overflow_scale.all_finite() && std::isinf(overflow_scale[0]), "scale can produce Inf");
    Vector overflow_axpy{std::numeric_limits<double>::max()};
    overflow_axpy.axpy(1.0, Vector{std::numeric_limits<double>::max()});
    expect_true(!overflow_axpy.all_finite() && std::isinf(overflow_axpy[0]), "axpy can produce Inf");

    Vector copy = listed;
    copy[0] = 9.0;
    expect_near(listed[0], 1.0, 1e-12, "copy independence");
    Vector moved = std::move(copy);
    expect_near(moved[0], 9.0, 1e-12, "move construction");

    return finish("MatCal::Linalg Vector tests");
}
