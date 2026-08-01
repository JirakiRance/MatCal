#include <limits>

#include "Matrix.hpp"
#include "test_support.hpp"

using namespace matcal_test;

namespace {

double residual_infinite_norm(const MatCal::Utils::Matrix& a,
                              const MatCal::Utils::Matrix& x,
                              const MatCal::Utils::Matrix& b) {
    double max_residual = 0.0;
    for (int i = 0; i < a.getRows(); ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < a.getCols(); ++j) {
            row_sum += a.get(i, j) * x.get(j, 0);
        }
        max_residual = std::max(max_residual, std::abs(row_sum - b.get(i, 0)));
    }
    return max_residual;
}

void expect_sor_converges(double omega) {
    using MatCal::Algorithm::Matrix::SOR;
    using MatCal::Utils::Matrix;

    Matrix a({{4.0, 1.0}, {2.0, 3.0}});
    Matrix b({{1.0}, {2.0}});
    auto result = SOR(a, b, omega, 1e-10, 500);
    expect_true(result.converged, "SOR converges for omega=" + std::to_string(omega));
    expect_near(result.root.get(0, 0), 0.1, 1e-7, "SOR solution x0 omega=" + std::to_string(omega));
    expect_near(result.root.get(1, 0), 0.6, 1e-7, "SOR solution x1 omega=" + std::to_string(omega));
    expect_true(residual_infinite_norm(a, result.root, b) < 1e-8, "SOR residual omega=" + std::to_string(omega));
}

} // namespace

int main() {
    using MatCal::Algorithm::Matrix::SOR;
    using MatCal::Utils::Matrix;

    expect_sor_converges(1.0);
    expect_sor_converges(1.2);
    expect_sor_converges(1.5);

    Matrix a({{4.0, 1.0}, {2.0, 3.0}});
    Matrix b({{1.0}, {2.0}});
    auto legacy_int = SOR(a, b, 1, 1e-10, 500);
    expect_true(legacy_int.converged, "legacy int SOR overload remains source-compatible");
    expect_true(residual_infinite_norm(a, legacy_int.root, b) < 1e-8, "legacy int SOR overload residual");

    expect_throw([&] { (void)SOR(a, b, 0.0); }, "SOR rejects omega=0");
    expect_throw([&] { (void)SOR(a, b, 2.0); }, "SOR rejects omega=2");
    expect_throw([&] { (void)SOR(a, b, std::numeric_limits<double>::quiet_NaN()); }, "SOR rejects NaN omega");
    expect_throw([&] { (void)SOR(a, b, std::numeric_limits<double>::infinity()); }, "SOR rejects Inf omega");

    Matrix divergent_a({{1.0, 3.0}, {3.0, 1.0}});
    Matrix divergent_b({{1.0}, {1.0}});
    auto non_converged = SOR(divergent_a, divergent_b, 1.0, 1e-14, 3);
    expect_true(!non_converged.converged, "SOR reports non-convergence when max iterations is too small");

    return finish("legacy SOR regression tests");
}
