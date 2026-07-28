int header_nonlinear_self_contained();
int header_least_squares_self_contained();
int header_polynomial_interpolation_self_contained();

int main() {
    return header_nonlinear_self_contained() +
           header_least_squares_self_contained() +
           header_polynomial_interpolation_self_contained();
}
