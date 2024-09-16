/**
 * @Author: Your name
 * @Date:   2024-09-15 11:09:40
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 14:55:16
 */
#include"head.h"
#include"basefun.h"
#include"DGSolver.h"
#pragma comment (lib, "libfftw3-3.lib")
#pragma comment (lib, "libfftw3f-3.lib")
#pragma comment (lib, "libfftw3l-3.lib")

int main()
{
	int Nx = 64;
	int Ny = 64;
	int Nz = 64;
	int dim, Dim;
	double c, dt, tau, gamma;
	fftw_plan plan_r2c_forward, plan_c2r_backward;

	int RealSize = Nx * Ny * Nz;
	int ImagSize = Nx * Ny * (Nz / 2 + 1);
	fftw_complex* Phi_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* Phi_hat_t = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_Nonlinear_terms_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_gradient = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_e2_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_e_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);

	double* Phi = fn_vec_init<double>(RealSize);
	double* Phi2 = fn_vec_init<double>(RealSize);
	double* Phi3 = fn_vec_init<double>(RealSize);
	double* Phi4 = fn_vec_init<double>(RealSize);
	double* LB_Nonlinear_terms = fn_vec_init<double>(RealSize);
	double* LB_e2 = fn_vec_init<double>(RealSize);
	double* LB_e = fn_vec_init<double>(RealSize);
	
	DGSolver solver_normal(Nx, Ny, Nz);
	solver_normal.set_para(&dim, &Dim, &c, &dt, &tau, &gamma);
	double** ProjMatrix = fn_mat_init<double>(dim, Dim);
	double** pk = fn_mat_init<double>(Nx * Ny * (Nz / 2 + 1), dim);
	double* pk2 = fn_vec_init<double>(Nx * Ny * (Nz / 2 + 1));
	double* kt = fn_vec_init<double>(Nx * Ny * (Nz / 2 + 1));
	double* ktmp = fn_vec_init<double>(Nx * Ny * (Nz / 2 + 1));
	double* kt_gradient = fn_vec_init<double>(Nx * Ny * (Nz / 2 + 1));

	solver_normal.ProjMatrix(ProjMatrix);
	solver_normal.get_kspace(ProjMatrix, pk, pk2, kt, ktmp, kt_gradient);
	solver_normal.fftw_plan_initial(&plan_r2c_forward, &plan_c2r_backward, Phi, Phi_hat);
	solver_normal.get_initial_value(Phi_hat);
	fftw_execute_dft_c2r(plan_c2r_backward, Phi_hat, Phi);

	solver_normal.get_initial_value(Phi_hat);// The previous step will change the value of Phi_hat
	
	multiply_xy(Phi, Phi, Phi2, Nx * Ny * Nz);
	multiply_xy(Phi, Phi2, Phi3, Nx * Ny * Nz);
	multiply_xy(Phi, Phi3, Phi4, Nx * Ny * Nz);
	
	for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
		LB_Nonlinear_terms[j1] = 0.5 * gamma * Phi2[j1] - Phi3[j1] / 6.0 - tau * Phi[j1];
	}
	
	fftw_execute_dft_r2c(plan_r2c_forward, LB_Nonlinear_terms, LB_Nonlinear_terms_hat);
	average_fourier(LB_Nonlinear_terms_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);
	double energy, gradient = 1.0, energy0 = 0.0, energy_err = 1.0;
	int iter = 0;
	double tol = 1e-14;
	int Iter_max = 3000;
	double* Gradient = fn_vec_init<double>(Iter_max);
	double* Energy = fn_vec_init<double>(Iter_max);
	while (energy_err > tol && iter < Iter_max)
	{
		iter = iter + 1;

		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			Phi_hat_t[j1][0] = (Phi_hat[j1][0] + dt * LB_Nonlinear_terms_hat[j1][0]) / ktmp[j1];
			Phi_hat_t[j1][1] = (Phi_hat[j1][1] + dt * LB_Nonlinear_terms_hat[j1][1]) / ktmp[j1];
		}
		Phi_hat_t[0][0] = 0.0;
		Phi_hat_t[0][1] = 0.0;
		
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			LB_gradient[j1][0] = fabs(Phi_hat_t[j1][0] - Phi_hat[j1][0]) / dt;
			LB_gradient[j1][1] = fabs(Phi_hat_t[j1][1] - Phi_hat[j1][1]) / dt;
		}
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			Phi_hat[j1][0] = Phi_hat_t[j1][0];
			Phi_hat[j1][1] = Phi_hat_t[j1][1];
		}

		gradient = Max(LB_gradient, Nx * Ny * (Nz / 2 + 1));
		Gradient[iter - 1] = gradient;
		
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			LB_e2_hat[j1][0] = kt[j1] * Phi_hat_t[j1][0];
			LB_e2_hat[j1][1] = kt[j1] * Phi_hat_t[j1][1];
		}
		fftw_execute_dft_c2r(plan_c2r_backward, Phi_hat_t, Phi);
		multiply_xy(Phi, Phi, Phi2, Nx * Ny * Nz);
		multiply_xy(Phi2, Phi, Phi3, Nx * Ny * Nz);
		multiply_xy(Phi3, Phi, Phi4, Nx * Ny * Nz);
		for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
			LB_Nonlinear_terms[j1] = 0.5 * gamma * Phi2[j1] - Phi3[j1] / 6.0 - tau * Phi[j1];
		}
		
		fftw_execute_dft_r2c(plan_r2c_forward, LB_Nonlinear_terms, LB_Nonlinear_terms_hat);
		average_fourier(LB_Nonlinear_terms_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);

		fftw_execute_dft_c2r(plan_c2r_backward, LB_e2_hat, LB_e2);
		multiply_xy(LB_e2, LB_e2, LB_e, Nx * Ny * Nz);
		for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
			LB_e[j1] = c * LB_e[j1] / 2.0 + tau * Phi2[j1] / 2.0 - gamma * Phi3[j1] / 6.0 + Phi4[j1] / 24.0;
		}
		fftw_execute_dft_r2c(plan_r2c_forward, LB_e, LB_e_hat);
		average_fourier(LB_e_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);

		energy = **LB_e_hat;
		*(Energy + iter - 1) = energy;
		energy_err = fabs(energy - energy0);
		energy0 = energy;
		std::cout << "iter = " << iter << " ";
		std::cout << "gradient = " << std::scientific << gradient << " ";
		std::cout << std::setprecision(15) << "energy = " << energy << " ";
		std::cout << std::setprecision(6) << "energy_err = " << std::scientific << energy_err << std::endl;
	}
	fn_mat_free<double>(ProjMatrix);
	fn_mat_free<double>(pk);
	fn_vec_free<double>(pk2);
	fn_vec_free<double>(kt);
	fn_vec_free<double>(ktmp);
	fn_vec_free<double>(kt_gradient);
	fn_vec_free<double>(Phi);
	fn_vec_free<double>(Phi2);
	fn_vec_free<double>(Phi3);
	fn_vec_free<double>(Phi4);
	fn_vec_free<double>(LB_Nonlinear_terms);
	fn_vec_free<double>(LB_e2);
	fn_vec_free<double>(LB_e);
	fn_vec_free<double>(Gradient);
	fn_vec_free<double>(Energy);
	fn_vec_free_cplx<fftw_complex>(Phi_hat);
	fn_vec_free_cplx<fftw_complex>(Phi_hat_t);
	fn_vec_free_cplx<fftw_complex>(LB_Nonlinear_terms_hat);
	fn_vec_free_cplx<fftw_complex>(LB_gradient);
	fn_vec_free_cplx<fftw_complex>(LB_e2_hat);
	fn_vec_free_cplx<fftw_complex>(LB_e_hat);
	return 0;
}