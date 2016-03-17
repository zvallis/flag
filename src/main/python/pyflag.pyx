
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

    ctypedef struct flagfb_parameters_t:
        int B_l
        int L
        int J_min_l
        int N
        int B_k
        int K
        int J_min_k
        int spin
        int upsample
        int reality
        int verbosity
        double knodes
        double k_interval

    ctypedef struct flagfb_wavscal_t:
        double complex *scal
        double complex *wav

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

    void flag_fourierbessel_mw_inverse(double complex *f, const double complex *flmn, const double *rnodes, const int Nrnodes, const double *knodes, const int Nknodes, const int L, const int spin, double complex *Flmr)

    void flag_fourierbessel_spherbessel_mapped_synthesis(double complex *f, const double complex *fn, const double *rnodes, int Nrnodes, const double *knodes, int Nknodes, int mapsize);

    int flag_k2kind(double k, double *knodes, double k_interval);

    double flag_kind2k(int k_ind, double *knodes, double k_interval);

    void fill_flaglet_parameters(flaglet_parameters_t *flaglet_parameters, const flagfb_parameters_t *parameters);

    double sjl(int l, double x);

    flagfb_wavscal_t flagfb_allocate_f_wav_scal(const flagfb_parameters_t *parameters)

    void flag_fbwavelet_analysis_lmnk(double complex *f_wav, double complex *f_scal, const double complex *flmp, const flagfb_parameters_t *parameters);

#----------------------------------------------------------------------------------------------------#

cdef extern from "stdlib.h":
    void free(void* ptr)

#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":
    int ssht_sampling_mw_n(int L);
    int ssht_sampling_mw_ntheta(int L);
    int ssht_sampling_mw_nphi(int L);

#----------------------------------------------------------------------------------------------------#

cdef extern from "flaglet.h":

    ctypedef struct flaglet_parameters_t:
        int B_l;
        int L;
        int J_min_l;
        int N;
        int B_p;
        int P;
        int J_min_p;
        int spin;
        int upsample;
        int reality;
        double tau;

    int flaglet_j_max(int L, int B);

    void flaglet_allocate_wav_lmp(double complex **wav_lmp, double **scal_lmp, const flaglet_parameters_t *parameters);

    void flaglet_wav_lmp(double complex *wav_lmp, double *scal_lmp, const flaglet_parameters_t *parameters);

    int flaglet_radial_bandlimit(int jp, const flaglet_parameters_t *parameters);

    int flaglet_angular_bandlimit(int jl, const flaglet_parameters_t *parameters);

#----------------------------------------------------------------------------------------------------#

def pyflag_k2kind(k, knodes, k_interval):
    kmin = min(knodes)
    kind = np.ceil((k-kmin)/k_interval)
    return int(kind)

def pyflag_kind2k(k_ind, knodes, k_interval):
    kmin = min(knodes)
    k = kmin + k_ind*k_interval
    return k

#----------------------------------------------------------------------------------------------------#

def healpy_lm(el, em, L):
    return em*(2*L-1-em)/2+el

def lmk2lmk_hp(np.ndarray[double complex, ndim=1, mode="c"] f_lmk not None, L, knodes, Nknodes, k_interval):
    kmapsize = L*L
    kmapsize_hp = L*(L+1)/2
    f_lmk_hp = np.zeros([Nknodes*L*(L+1)/2,],dtype=complex)
    for k_ind from 0<= k_ind <Nknodes:
        k = pyflag_kind2k(k_ind,knodes,k_interval)
        for el from 0 <= el < L:
            for em from 0 <= em <= el:
                f_lmk_hp[healpy_lm(el, em, L) + k_ind*kmapsize_hp] = f_lmk[ el * el + el + em + k_ind*kmapsize]
    return f_lmk_hp

def lmk_hp2lmk(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, L, knodes, Nknodes, k_interval):
    kmapsize = L*L
    kmapsize_hp = L*(L+1)/2
    f_lm = np.zeros([Nknodes*L*L,], dtype=complex)
    for k_ind from 0<= k_ind <Nknodes:
        k = pyflag_kind2k(k_ind,knodes,k_interval)
        for el from 0 <= el < L:
            for em from 0 <= em <= el:
                f_lm[ el * el + el - em + k_ind*kmapsize] = pow(-1.0, -em) * (flm_hp[healpy_lm(el, em, L) + k_ind*kmapsize_hp] ).conjugate()
                f_lm[ el * el + el + em + k_ind*kmapsize] = flm_hp[healpy_lm(el, em, L) + k_ind*kmapsize_hp]
    return f_lm

#----------------------------------------------------------------------------------------------------#

def pyflag_fourierbessel_mw_inverse(np.ndarray[double complex, ndim=1, mode="c"] f, np.ndarray[double complex, ndim=1, mode="c"] flmn, np.ndarray[double, ndim=1, mode="c"] rnodes, Nrnodes, np.ndarray[double, ndim=1, mode="c"] knodes, Nknodes, L, spin, np.ndarray[double complex, ndim=1, mode="c"] Flmr):
    print("pyflag: Nrnodes = " + str(Nrnodes) + " Nknodes = " + str(Nknodes) + " L = " + str(L) + " spin = " + str(spin))
    flag_fourierbessel_mw_inverse(<double complex*> np.PyArray_DATA(f), <double complex*> np.PyArray_DATA(flmn), <double*> np.PyArray_DATA(rnodes), Nrnodes, <double*> np.PyArray_DATA(knodes), Nknodes, L, spin, <double complex*> np.PyArray_DATA(Flmr))
    if Flmr is None:
        print("\nFlmr is None \n")

#----------------------------------------------------------------------------------------------------#

#def flag_n_scal(const flagfb_parameters_t *parameters):
#    J_min = parameters.J_min_l
#    L = parameters.L
#    if (parameters.upsample):
#        bandlimit = parameters.L
#    else:
#        bandlimit = min(s2let_bandlimit(J_min-1, parameters), L)
#    #Currently only for MW sampling
#    ntheta = ssht_sampling_mw_ntheta(L)
#    nphi = ssht_sampling_mw_nphi(L)
#    return nphi * ntheta;

#def flag_n_wav(const flagfb_parameters_t *parameters):
#    so3_parameters_t so3_parameters = {}
#    fill_so3_parameters(&so3_parameters, parameters)
#
#    s2let_parameters_t s2let_parameters = {}
#    s2let_parameters.B = parameters.B

#    L = parameters.L;
#    J_min = parameters.J_min_l;
#    J = s2let_j_max(parameters);
#    bandlimit = L;
#    total = 0;
#    for j in range (J_min,J):
#        if (not parameters.upsample):
#            bandlimit = min(s2let_bandlimit(j, parameters), L)
#            so3_parameters.L = bandlimit
#        total = total + so3_sampling_f_size(&so3_parameters)
#    print("total")
#    print(total)
#    return total

#----------------------------------------------------------------------------------------------------#

# Currently only for MW, no MW_SS implementation
def pyso3_sampling_f_size(L,N):
    nalpha = 2*L-1
    nbeta = L
    # Assuming steerable = 1
    ngamma = N
    return nalpha * nbeta * ngamma

def pyflag_define_f_wav(int L, double K, int N, int B_l, double B_k):
    L = parameters.L
    K = parameters.K
    N = parameters.N
    J_l = flaglet_j_max(L, parameters.B_l);
    J_k = flaglet_j_max(K, parameters.B_k);
    jp, jl, total = 0;
    Nj = min(N,L);
    bandlimit_k = K;
    bandlimit_l = L;
    cdef flaglet_parameters_t flaglet_parameters = {}
    fill_flaglet_parameters(&flaglet_parameters, parameters)
    for jk in range (parameters.J_min_k,J_k+1):
        if (not parameters.upsample):
            bandlimit_k = min(flaglet_radial_bandlimit(jk, &flaglet_parameters), K)
        for jl in range (parameters.J_min_l,J_l+1):
            if (not parameters.upsample):
                bandlimit_l = min(flaglet_angular_bandlimit(jl, &flaglet_parameters), L)
                so3_L = bandlimit_l
                Nj = min(N,bandlimit_l)
                Nj = Nj (Nj+N)%2
                so3_N = Nj
            L0 = np.ceil(pow(parameters.B_l, jl-1))
            total = total + pyso3_sampling_f_size(so3_L,so3_N) * bandlimit_k
    #*f_wav = (complex double*)calloc( total, sizeof(complex double));
    #*f_scal = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
    f_scal = np.array((L,2*L-1,K),'complex')
    f_wav = np.array((N,L,2*L-1,bandlimit_k),'complex')
    return [f_scal,f_wav]

#----------------------------------------------------------------------------------------------------#

def pyflag_fbanalysis_lmk2wav(np.ndarray[double complex, ndim=1, mode="c"] flmk_hp not None, B_l, B_k, L, J_min, N, spin, upsample, knodes, k_max, k_interval):

    cdef flagfb_parameters_t parameters = {};
    parameters.B_l = B_l;
    parameters.B_k = B_k
    parameters.L = L;
    parameters.J_min_l = J_min;
    parameters.N = N;
    parameters.K = k_max
    parameters.spin = spin;
    parameters.upsample = upsample;
    parameters.reality = 0
    parameters.verbosity = 0
    parameters.knodes = knodes
    parameters.k_interval = k_interval
    print(knodes)

    #going to try moving to c instead
    #f_scal = np.zeros([pyflag_n_scal(&parameters),], dtype=complex)
    #f_wav = np.zeros([pyflag_n_wav(&parameters),], dtype=complex)
    #f_wavscal = flagfb_allocate_f_wav_scal(&parameters)
    #f_scal = f_wavscal.scal
    #f_wav = f_wavscal.wav
    #flagfb_allocate_f_wav

    [f_scal,f_wav] = pyflag_define_f_wav(parameters)
    #f_lmp = lmk_hp2lmk(flmk_hp, L, knodes, Nknodes, k_interval)
    print("f_scal shape ")
    print(f_scal.shape)
    print("f_wav shape ")
    print(f_wav.shape)

    #s2let_analysis_lm2wav(<double complex*> np.PyArray_DATA(f_wav), <double complex*> np.PyArray_DATA(f_scal), <const double complex*> np.PyArray_DATA(f_lm), &parameters);

    #flag_fbwavelet_analysis_lmnk(<double complex*> np.PyArray_DATA(f_wav), <double complex*> np.PyArray_DATA(f_scal), <const double complex*> np.PyArray_DATA(f_lmp), &parameters)

    #return f_wav, f_scal
