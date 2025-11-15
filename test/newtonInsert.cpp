#include <iostream>
#include <vector>
#include "..\src\NI.hpp"   // Your class
using namespace std;
using namespace MatCal::Utils;

// A test function to validate interpolation correctness
double f(double x) {
    return x*x*x + 2*x - 5;  // You may change to any function you like
}

int main() {

    //------------------------------------------------------------
    // 1. Construct equally spaced sample points
    //------------------------------------------------------------
    double x0 = 0.0;
    double h  = 1.0;

    vector<pair<double,double>> data;
    for(int i=0; i<3; i++){          // Generate x = 0,1,2
        double x = x0 + i*h;
        data.push_back({x, f(x)});
    }

    //------------------------------------------------------------
    // 2. Build the interpolator
    //------------------------------------------------------------
    NewtonInsert_Finite interp(h, data);

    cout << "=== Initial Interpolating Polynomial: degree = " 
         << interp.getDegree() << " ===\n";
    
    for(double X=0; X<=2; X+=0.5){
        cout << "P(" << X << ") = " << interp.calculate(X) 
             << "\t   f(" << X << ") = " << f(X) << endl;
    }

    //------------------------------------------------------------
    // 3. Print the finite difference table
    //------------------------------------------------------------
    cout << "\n=== Finite Difference Table ===\n";
    const auto& sheet = interp.getSheet();

    for(int i=0; i<=interp.getDegree(); i++){
        for(int j=0; j<=interp.getDegree(); j++){
            double v = sheet.get(i, j);
            if(j > i) v = 0;      // Only lower triangular part is valid
            cout << v << "\t";
        }
        cout << endl;
    }

    //------------------------------------------------------------
    // 4. Insert a new point (x = 3)
    //------------------------------------------------------------
    double x_new = 3;
    double y_new = f(x_new);

    cout << "\nInsert new point: x = " << x_new 
         << ", y = " << y_new << endl;

    interp.insertNewTerm(y_new);

    cout << "=== Polynomial After Insertion: degree = " 
         << interp.getDegree() << " ===\n";

    for(double X=0; X<=3; X+=0.5){
        cout << "P(" << X << ") = " << interp.calculate(X) 
             << "\t   f(" << X << ") = " << f(X) << endl;
    }

    //------------------------------------------------------------
    // 5. Check interpolation error
    //------------------------------------------------------------
    cout << "\n=== Error Check ===\n";
    for(double X=0; X<=3; X+=0.5){
        double err = interp.calculate(X) - f(X);
        cout << "X = " << X << ", error = " << err << endl;
    }

    return 0;
}
