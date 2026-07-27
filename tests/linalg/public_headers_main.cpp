int matcal_linalg_header_vector_probe();
int matcal_linalg_header_dense_matrix_probe();
int matcal_linalg_header_solver_types_probe();
int matcal_linalg_header_dense_solver_probe();

int main() {
    return matcal_linalg_header_vector_probe() +
           matcal_linalg_header_dense_matrix_probe() +
           matcal_linalg_header_solver_types_probe() +
           matcal_linalg_header_dense_solver_probe();
}
