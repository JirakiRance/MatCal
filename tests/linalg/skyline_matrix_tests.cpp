#include <limits>

#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

int main() {
    using MatCal::Linalg::SymmetricSkylineMatrix;
    using MatCal::Linalg::Vector;

    SymmetricSkylineMatrix empty;
    expect_true(empty.size() == 0 && empty.storage_size() == 0, "default skyline is empty");

    SymmetricSkylineMatrix diagonal(3);
    expect_true(diagonal.storage_size() == 3, "diagonal skyline storage");
    diagonal.set(0, 0, 1.0);
    diagonal.add(1, 1, 2.0);
    diagonal.set(2, 2, 3.0);
    expect_near(diagonal.get(1, 1), 2.0, 1e-12, "diagonal add");
    expect_true(!diagonal.stores(2, 0), "diagonal profile excludes offdiagonal");
    expect_near(diagonal.get(2, 0), 0.0, 1e-12, "profile outside get returns structural zero");
    expect_throw([&] { diagonal.set(2, 0, 1.0); }, "set outside profile throws");
    expect_throw([&] { diagonal.add(0, 2, 1.0); }, "add outside symmetric profile throws");
    expect_throw([&] { (void)diagonal.get(3, 0); }, "get out of range throws");

    SymmetricSkylineMatrix skyline = SymmetricSkylineMatrix::from_profile({0, 0, 0});
    expect_true(skyline.storage_size() == 6, "full 3x3 lower profile storage");
    skyline.set(0, 0, 4.0);
    skyline.set(1, 0, 1.0);
    skyline.set(1, 1, 3.0);
    skyline.set(2, 0, 0.0);
    skyline.set(2, 1, 2.0);
    skyline.set(2, 2, 5.0);
    expect_near(skyline.get(0, 1), 1.0, 1e-12, "symmetric get");
    expect_near(skyline.get(1, 2), 2.0, 1e-12, "symmetric get upper");
    expect_near(skyline.matrix_scale(), 5.0, 1e-12, "matrix scale");
    expect_true(skyline.all_finite(), "finite skyline");

    Vector product = skyline.multiply(Vector{1.0, 2.0, 3.0});
    expect_near(product[0], 6.0, 1e-12, "skyline matvec row0");
    expect_near(product[1], 13.0, 1e-12, "skyline matvec row1");
    expect_near(product[2], 19.0, 1e-12, "skyline matvec row2");
    expect_throw([&] { (void)skyline.multiply(Vector{1.0}); }, "matvec dimension mismatch");

    auto from_pairs = SymmetricSkylineMatrix::from_symmetric_positions(4, {{0, 0}, {3, 1}, {2, 2}});
    expect_true(from_pairs.stores(3, 1), "pair profile stores requested position");
    expect_true(from_pairs.stores(3, 2), "pair profile stores skyline gap");
    expect_true(!from_pairs.stores(2, 0), "unrelated profile remains sparse");
    expect_true(from_pairs.storage_size() == 6, "pair profile storage includes skyline gap");

    expect_throw([] { SymmetricSkylineMatrix::from_profile({1}); }, "invalid first column rejected");
    expect_throw([] { (void)SymmetricSkylineMatrix::from_symmetric_positions(2, {{0, 2}}); }, "pair out of range rejected");

    SymmetricSkylineMatrix non_finite(1);
    non_finite.set(0, 0, std::numeric_limits<double>::infinity());
    expect_true(!non_finite.all_finite(), "finite check detects Inf");

    return finish("MatCal::Linalg SymmetricSkylineMatrix tests");
}
