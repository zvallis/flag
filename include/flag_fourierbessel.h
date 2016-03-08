#ifndef FLAG_FOURIERBESSEL
#define FLAG_FOURIERBESSEL

#include <complex.h>

void flag_fourierbessel_mw_inverse(complex double *f,
		const complex double *flmn,
		const double *nodes, int Nnodes,
		int L, double tau, int N, int spin);

void flag_fourierbessel_spherbessel_mapped_synthesis(complex double *f, const complex double *fn, const double *nodes, int Nnodes, double tau, int K, int mapsize);

int flag_k2kind(double k, double *knodes, double k_interval);

double flag_kind2k(int k_ind, double *knodes, double k_interval);

double sjl(int l, double x);

#endif