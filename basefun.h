#pragma once
#include"head.h"

template <typename T>
T* fn_vec_init_cplx(int n)
{
	T* src = (T*)fftw_malloc(sizeof(T) * n);
	return src;
}

template <typename T>
T* fn_vec_init(int n)
{
	T* src = (T*)malloc(sizeof(T) * n);
	return src;
}

template <typename T>
T** fn_mat_init_cplx(int m, int n)
{
	T** src = (T**)fftw_malloc(sizeof(T*) * m);
	T* dataBlock = (T*)fftw_malloc(m * n * sizeof(T));
	for (size_t i = 0; i < m; ++i) {
		src[i] = dataBlock + i * n;
	}
	return src;
}

template <typename T>
T** fn_mat_init(int m, int n)
{
	T** src = (T**)malloc(sizeof(T*) * m);
	T* dataBlock = (T*)malloc(m * n * sizeof(T));
	for (size_t i = 0; i < m; ++i) {
		src[i] = dataBlock + i * n;
	}
	return src;
}

template <typename T>
void fn_vec_free_cplx(T* src)
{
	fftw_free(src);
}

template <typename T>
void fn_vec_free(T* src)
{
	free(src);
}

template <typename T>
void fn_mat_free_cplx(T** src)
{
	fftw_free(src[0]);
	fftw_free(src);
}

template <typename T>
void fn_mat_free(T** src)
{
	free(src[0]);
	free(src);
}

inline void multiply_xy(double* in_x, double* in_y, double* out, int length)
{
	for (int j = 0; j < length; j++) {
		out[j] = in_x[j] * in_y[j];
	}
}

inline void average_fourier(fftw_complex* u, int length_u, int average_M)
{
	for (int j = 0; j < length_u; j++) {
		u[j][0] = u[j][0] / (average_M * 1.0);
		u[j][1] = u[j][1] / (average_M * 1.0);
	}
}

inline double Max(fftw_complex* array, int length)
{
	double max = 0.0;
	for (int j = 0; j < length; j++) {
		double temp = sqrt(array[j][0] * array[j][0] + array[j][1] * array[j][1]);
		if (max < temp) max = temp;
	}
	return max;
}