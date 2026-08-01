int matcal_linalg_header_vector_probe();
int matcal_linalg_header_dense_matrix_probe();
int matcal_linalg_header_solver_types_probe();
int matcal_linalg_header_dense_solver_probe();
int header_iterative_solvers_self_contained();
int header_eigen_solvers_self_contained();
int matcal_linalg_header_symmetric_skyline_probe();
int matcal_linalg_header_skyline_ldlt_probe();

int main() {
    return matcal_linalg_header_vector_probe() +
           matcal_linalg_header_dense_matrix_probe() +
           matcal_linalg_header_solver_types_probe() +
           matcal_linalg_header_dense_solver_probe() +
           header_iterative_solvers_self_contained() +
           header_eigen_solvers_self_contained() +
           matcal_linalg_header_symmetric_skyline_probe() +
           matcal_linalg_header_skyline_ldlt_probe();
}
