int matcal_roots_header_self_contained();
int matcal_linear_interpolator_header_self_contained();
int matcal_cubic_spline_header_self_contained();

int main() {
    return matcal_roots_header_self_contained() +
           matcal_linear_interpolator_header_self_contained() +
           matcal_cubic_spline_header_self_contained();
}
