#include <iostream>
#include <vector>
#include <cmath>
#include "..\src\Insert.hpp"

using namespace MatCal::Algorithm::Insert;

void testLagrange() {
    std::cout << "=== LagrangeInsert Test ===\n";
    std::vector<std::pair<double,double>> data = {
        {0, 1}, {1, 2}, {2, 5}, {3, 10}
    };
    LagrangeInsert lag(data);

    std::cout << "Degree: " << lag.getDegree() << "\n";
    for (double x = 0; x <= 3; x += 0.5) {
        std::cout << "x=" << x << ", y=" << lag.calculate(x) << "\n";
    }
    // Reconstruct test
    std::vector<std::pair<double,double>> newData = {
        {0, 0}, {1, 1}, {2, 8}
    };
    lag.reconstruct(newData);
    std::cout << "After reconstruct:\n";
    for (double x = 0; x <= 2; x += 0.5)
        std::cout << "x=" << x << ", y=" << lag.calculate(x) << "\n";
}

void testNewtonQuotient() {
    std::cout << "\n=== NewtonInsert_Quotient Test ===\n";
    std::vector<std::pair<double,double>> data = {
        {0,1},{1,2},{2,5},{3,10}
    };
    NewtonInsert_Quotient nq(data);

    std::cout << "Degree: " << nq.getDegree() << "\n";
    for (double x = 0; x <= 3; x += 0.5)
        std::cout << "x=" << x << ", y=" << nq.calculate(x) << "\n";

    // Test insertNewTerm
    try {
        nq.insertNewTerm(4, 17);
        std::cout << "After insertNewTerm x=4, y=17:\n";
        std::cout << "x=4, y=" << nq.calculate(4) << "\n";
    } catch(const std::exception &e) {
        std::cout << "Exception: " << e.what() << "\n";
    }

    // Test out-of-range insert
    try {
        NewtonInsert_Quotient nq_empty;
        nq_empty.insertNewTerm(1,1);
    } catch(const std::exception &e) {
        std::cout << "Expected exception: " << e.what() << "\n";
    }
}

void testNewtonFinite() {
    std::cout << "\n=== NewtonInsert_Finite Test ===\n";
    std::vector<double> yData = {1,2,5,10};
    double h = 1.0;
    double x0 = 0.0;
    NewtonInsert_Finite nf(h, x0, yData);

    std::cout << "Degree: " << nf.getDegree() << "\n";
    for (double x = 0; x <= 3; x += 0.5)
        std::cout << "x=" << x << ", y=" << nf.calculate(x) << "\n";

    // Test insertNewTerm
    try {
        nf.insertNewTerm(17);
        std::cout << "After insertNewTerm Y=17:\n";
        std::cout << "x=4, y=" << nf.calculate(4) << "\n";
    } catch(const std::exception &e) {
        std::cout << "Exception: " << e.what() << "\n";
    }

    // Test out-of-range insert
    try {
        std::vector<double> tiny = {1};
        NewtonInsert_Finite nf2(h,0,tiny);
        nf2.insertNewTerm(2);
    } catch(const std::exception &e) {
        std::cout << "Expected exception: " << e.what() << "\n";
    }
}

int main() {
    testLagrange();
    testNewtonQuotient();
    testNewtonFinite();
    return 0;
}
