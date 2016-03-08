
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "smk_shbessel.h":
    double sjl(int l, double x);

#----------------------------------------------------------------------------------------------------#

cdef extern from "flag.h":

    int ssht_fr_size_mw(int L);

    int ssht_flm_size(int L);

    int flag_core_flmn_size(int L, int N);

    int flag_core_f_size_mw(int L, int N);

    void flag_core_allocate_flmn(double complex **flmn, int L, int N);

    void flag_core_allocate_f(double complex **f, int L, int N);

    void flag_core_allocate_f_real(double **f, int L, int N);

    void flag_core_analysis(double complex *flmn,
        const double complex *f,
        int L, double tau, int N, int spin);

    void flag_core_synthesis(double complex *f,
        const double complex *flmn,
        const double *nodes, int Nnodes,
        int L, double tau, int N, int spin);

    void flag_core_analysis_real(double complex *flmn,
        const double *f, int L, double tau, int N);

    void flag_core_synthesis_real(double *f,
        const double complex *flmn,
        const double *nodes, int Nnodes,
        int L, double tau, int N);

    double j_ell(double X, int l);

    void flag_spherbessel_basis(double *jell, const int ell, const double *nodes, int Nnodes);

    void flag_core_fourierbessel_analysis(double complex *flmn,
        const double complex *f,
        int L, double tau, int N);

    void flag_core_fourierbessel_synthesis(double complex *f,
        const double complex *flmn,
        const double *nodes, int Nnodes,
        int L, double tau, int N);


    void allocate_ssht_sampling(double **thetas, double **phis, int L);

    void flag_sampling_allocate(double **rs, double **thetas, double **phis, double **laguweights, int L, int N);

    void ssht_sampling_mw(double *thetas, double *phis, int L);

    void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double tau, int L, int N);

    void flag_spherlaguerre2spherbessel(double *flk, const double *fn, double *kvalues, int Nk, int N, int ell, double tau);

    void flag_fourierlaguerre2fourierbessel(double complex *flmk, double complex *flmn, double *kvalues, int Nk, int N, int L, double tau);

    void flag_sbesselslag(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

    void flag_mulk(double *mulk, int n, int ell, double k, double tau);

    void bubbleSort(double numbers[], int array_size);

    void flag_sbesselslag_backup(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

    void flag_sbesselslag_backup2(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

    double eval_laguerre(double z, int n, int alpha);

    double eval_laguerre_rescaled(double z, int n, int alpha, double normfac);

    long factorial(int n);

    long factorial_range(int min, int max);

    void flag_spherlaguerre_quadrature(double *roots, double *weights, int N, int alpha);

    double flag_spherlaguerre_tau(double R, int N);

    double flag_spherlaguerre_Rmax(int N);

    void flag_spherlaguerre_sampling(double *nodes, double *weights, double tau, int N);

    void flag_spherlaguerre_allocate_sampling(double **nodes, double **weights, int N);

    void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, double tau, int N);

    void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N);

    void flag_spherlaguerre_synthesis_gen(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N, int alpha);

    void flag_spherlaguerre_mapped_analysis(double complex *fn, const double complex *f, const double *nodes, const double *weights, double tau, int N, int mapsize);

    void flag_spherlaguerre_mapped_synthesis(double complex *f, const double complex *fn, const double *nodes, int Nnodes, double tau, int N, int mapsize);

    void flag_spherlaguerre_basis(double *KN, const int N, const double *nodes, int Nnodes, double tau);

    void flag_spherbessel_sampling(double *nodes, double *weights, double R, int N);

    void flag_sbesselslag_backup3(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

    void flag_fourierbessel_mw_inverse(double complex *f, const double complex *flmn, const double *rnodes, const double *knodes, const int L, const int spin, const int Nrnodes, const int Nknodes)

    void flag_fourierbessel_spherbessel_mapped_synthesis(double complex *f, const double complex *fn, const double *rnodes, int Nrnodes, const double *knodes, int Nknodes, int mapsize);

    int flag_k2kind(double k, double *knodes, double k_interval);

    double flag_kind2k(int k_ind, double *knodes, double k_interval);

    double sjl(int l, double x);

    #void test(const double *arr, const int *L, const int *spin, const int *Nrnodes, const int *Nknodes);

#----------------------------------------------------------------------------------------------------#

cdef extern from "stdlib.h":
    void free(void* ptr)

#----------------------------------------------------------------------------------------------------#

def pyflag_k2kind(k, knodes, k_interval):
    kmin = min(knodes)
    kind = (k-kmin)/k_interval
    return int(kind)

def pyflag_kind2k(k_ind, knodes, k_interval):
    kmin = min(knodes)
    k = kmin + k_ind*k_interval
    return k

#----------------------------------------------------------------------------------------------------#

def pyflag_fourierbessel_mw_inverse(np.ndarray[double complex, ndim=1, mode="c"] f, np.ndarray[double complex, ndim=1, mode="c"] flmn, np.ndarray[double, ndim=1, mode="c"] rnodes, Nrnodes, np.ndarray[double, ndim=1, mode="c"] knodes, Nknodes, L, spin):
    print("pyflag: Nrnodes = " + str(Nrnodes) + " Nknodes = " + str(Nknodes) + " L = " + str(L) + " spin = " + str(spin))
    flag_fourierbessel_mw_inverse(<double complex*> np.PyArray_DATA(f), <double complex*> np.PyArray_DATA(flmn), <double*> np.PyArray_DATA(rnodes), <double*> np.PyArray_DATA(knodes), L, spin, Nrnodes, Nknodes)


def test_pass(np.ndarray[double, ndim=1, mode="c"] arr, Nrnodes, Nknodes, L, spin):
    print("TEST pyflag: L = " + str(L) + " Nrnodes = " + str(Nrnodes) + " Nknodes = " + str(Nknodes) + " spin = " + str(spin))
    print(arr)
    arrL = np.array([L])
    arrspin = np.array([spin])
    arrNrnodes = np.array([Nrnodes])
    arrNknodes = np.array([Nknodes])
    #test( <double*> np.PyArray_DATA(arr), <int*> np.PyArray_DATA(arrL), <int*> np.PyArray_DATA(arrspin), <int*> np.PyArray_DATA(arrNrnodes), <int*> np.PyArray_DATA(arrNknodes))
