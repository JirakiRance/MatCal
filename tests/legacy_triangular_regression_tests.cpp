#include "test_support.hpp"

using namespace matcal_test;

namespace {

void upper_copy_and_solve_are_safe() {
    using MatCal::Utils::Matrix;
    using MatCal::Utils::UpperTriangularMatrix;

    UpperTriangularMatrix upper(3);
    upper.set(0, 0, 2.0);
    upper.set(0, 1, -1.0);
    upper.set(0, 2, 1.0);
    upper.set(1, 1, 3.0);
    upper.set(1, 2, 2.0);
    upper.set(2, 2, 4.0);

    expect_true(upper.getRows() == 3 && upper.getCols() == 3, "upper original dimensions");
    expect_true(upper.getData().size() == 6, "upper original storage size");

    UpperTriangularMatrix copy(upper);
    expect_true(copy.getRows() == 3 && copy.getCols() == 3, "upper non-const copy dimensions");
    expect_true(copy.getData().size() == 6, "upper non-const copy storage size");
    expect_near(copy.get(0, 2), 1.0, 1e-12, "upper non-const copy get");
    copy.set(0, 2, 5.0);
    expect_near(copy.get(0, 2), 5.0, 1e-12, "upper non-const copy set");
    expect_near(upper.get(0, 2), 1.0, 1e-12, "upper copy is independent");

    const UpperTriangularMatrix const_upper(upper);
    UpperTriangularMatrix const_copy(const_upper);
    expect_true(const_copy.getRows() == 3 && const_copy.getCols() == 3, "upper const copy dimensions");
    expect_true(const_copy.getData().size() == 6, "upper const copy storage size");

    UpperTriangularMatrix assigned;
    assigned = upper;
    expect_true(assigned.getRows() == 3 && assigned.getCols() == 3, "upper assignment dimensions");
    expect_near(assigned.get(1, 2), 2.0, 1e-12, "upper assignment data");

    auto lower = upper.transpose();
    expect_true(lower->getRows() == 3 && lower->getCols() == 3, "upper transpose dimensions");
    expect_near(lower->get(2, 0), 1.0, 1e-12, "upper transpose data");

    Matrix rhs({{3.0}, {12.0}, {12.0}});
    Matrix x = as_matrix(upper.solve(rhs));
    expect_near(x.get(0, 0), 1.0, 1e-12, "upper solve x0");
    expect_near(x.get(1, 0), 2.0, 1e-12, "upper solve x1");
    expect_near(x.get(2, 0), 3.0, 1e-12, "upper solve x2");
}

void lower_copy_and_solve_are_safe() {
    using MatCal::Utils::LowerTriangularMatrix;
    using MatCal::Utils::Matrix;

    LowerTriangularMatrix lower(3);
    lower.set(0, 0, 2.0);
    lower.set(1, 0, -1.0);
    lower.set(1, 1, 3.0);
    lower.set(2, 0, 1.0);
    lower.set(2, 1, 2.0);
    lower.set(2, 2, 4.0);

    expect_true(lower.getRows() == 3 && lower.getCols() == 3, "lower original dimensions");
    expect_true(lower.getData().size() == 6, "lower original storage size");

    LowerTriangularMatrix copy(lower);
    expect_true(copy.getRows() == 3 && copy.getCols() == 3, "lower non-const copy dimensions");
    expect_true(copy.getData().size() == 6, "lower non-const copy storage size");
    expect_near(copy.get(2, 0), 1.0, 1e-12, "lower non-const copy get");
    copy.set(2, 0, 5.0);
    expect_near(copy.get(2, 0), 5.0, 1e-12, "lower non-const copy set");
    expect_near(lower.get(2, 0), 1.0, 1e-12, "lower copy is independent");

    const LowerTriangularMatrix const_lower(lower);
    LowerTriangularMatrix const_copy(const_lower);
    expect_true(const_copy.getRows() == 3 && const_copy.getCols() == 3, "lower const copy dimensions");
    expect_true(const_copy.getData().size() == 6, "lower const copy storage size");

    LowerTriangularMatrix assigned;
    assigned = lower;
    expect_true(assigned.getRows() == 3 && assigned.getCols() == 3, "lower assignment dimensions");
    expect_near(assigned.get(2, 1), 2.0, 1e-12, "lower assignment data");

    auto upper = lower.transpose();
    expect_true(upper->getRows() == 3 && upper->getCols() == 3, "lower transpose dimensions");
    expect_near(upper->get(0, 2), 1.0, 1e-12, "lower transpose data");

    Matrix rhs({{2.0}, {5.0}, {17.0}});
    Matrix x = as_matrix(lower.solve(rhs));
    expect_near(x.get(0, 0), 1.0, 1e-12, "lower solve x0");
    expect_near(x.get(1, 0), 2.0, 1e-12, "lower solve x1");
    expect_near(x.get(2, 0), 3.0, 1e-12, "lower solve x2");
}

} // namespace

int main() {
    upper_copy_and_solve_are_safe();
    lower_copy_and_solve_are_safe();
    return finish("legacy triangular regression tests");
}
