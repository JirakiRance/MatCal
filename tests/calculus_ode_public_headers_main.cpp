int matcal_calculus_header_self_contained();
int matcal_ode_header_self_contained();

int main() {
    return matcal_calculus_header_self_contained() + matcal_ode_header_self_contained();
}
