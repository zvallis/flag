//Changes by Zoe Vallis, based on code by Boris Leistedt and Jason McEwen

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>
//#include <assert.h>

#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"
#include "ssht_sampling.h"
#include "ssht_core.h"
#include "smk_shbessel.h"
#include "flag_core.h"

/*void test(const double *arr, const int *L, const int *spin, const int *Nrnodes, const int *Nknodes)
{
    printf("TEST inside flag_fb: L = %i, Nrnodes %i, Nknodes %i, spin %i\n",L[0],Nrnodes[0],Nknodes[0],spin[0]);
    printf("arr = %f %f %f %f \n", arr[0],arr[1],arr[2],arr[3]);
    return 0;
}*/

// Based on flag_core_synthesis, computes f from a given flmn
void flag_fourierbessel_mw_inverse(complex double *f, const complex double *flmn, const double *rnodes, const double *knodes, const int L, const int spin, const int Nrnodes, const int Nknodes)
{
    //assert(L > 0);
	//assert(N > 1);
	//assert(nodes != NULL);
	//const int alpha = ALPHA;

    printf("FLAG_FB start: Nrnodes %i, Nknodes %i, L %i, spin %i\n",Nrnodes,Nknodes,L,spin);
	int verbosity = 1; // Level of how much detail to print about ongoing process
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *Flmr;
	flag_core_allocate_flmn(&Flmr, L, Nrnodes);
	printf("> Mapped spherical Laguerre transform...");fflush(NULL);
    //Using variation that integrates over k, as opposed to Fourier-Laguerre which sums over p
    flag_fourierbessel_spherbessel_mapped_synthesis(Flmr, flmn, rnodes, Nrnodes, knodes, Nknodes, flmsize);
	printf("done\n");
    printf("FLAG_FB outside loop L = %i\n",L);

	for (n = 0; n < Nrnodes; n++){
        printf("FLAG_FB inside loop L = %i\n",L);
		printf("> Synthesis: layer %i on %i\n",n+1,Nrnodes);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, Flmr + offset_lm, L, spin, dl_method, verbosity);
	}

    free(Flmr);

}

// Based on flag_spherlaguerre_mapped_synthesis, changed to integrate over k with additional Bessel function, rather than summing over p
void flag_fourierbessel_spherbessel_mapped_synthesis(complex double *f, const complex double *fn, const double *rnodes, int Nrnodes, const double *knodes, int Nknodes, int mapsize)
{
	//assert(Nrnodes > 1);
    //assert(Nknodes > 1);
	//assert(mapsize > 1);
    int i, l, x, offset_n, offset_i;
    double r, k, k_min = knodes[0], k_max = knodes[Nknodes-1], k_interval = (k_max-k_min)/(Nknodes-1), s, S = 30; // S is the number of strips to divide the integration into;
    double s_width = ((k_max-k_min)/S)/6.0;
	double *bessel_integral = (double*)calloc(Nknodes, sizeof(double));

	for (i = 0; i < Nrnodes; i++)
	{
		r = rnodes[i];
        // Integrate over k
		for (k=k_min; k<k_max; k+=k_interval)
		{
            int k_ind = flag_k2kind(k,knodes,k_interval);
            offset_n = k_ind * mapsize;
                for(l=0; l<mapsize; l++)
                {
                    // Weddle's rule for integration
                    for (s=1; s<=S; s++)
                    {
                        int c = 0;
                        double h = s_width/6.0, y[7] = { 0 };
                        for (x=k; x<=k+s_width; x+=h)
                        {
                            y[c] = pow(x,2) * sjl(l,k*r) * fn[l+offset_n];
                            c += 1;
                        }
                        *bessel_integral += 3.0/10.0 * h * (y[0] + 5 * y[1] + y[2] + 6 * y[3] + y[4] + 5 * y[5] + y[6]);
                    }
                }
		}
        offset_i = i * mapsize;
        for(l=0; l<mapsize; l++)
			{
                f[l+offset_i] = *bessel_integral;
            }
    }

	free(bessel_integral);

}

//Convert from k value to corresponding k index
int flag_k2kind(double k, double *knodes, double k_interval){
    double kmin = knodes[0];
    int kind = (int)((k-kmin)/k_interval);
    return kind;
}

// Convert from k index to corresponding k value
double flag_kind2k(int k_ind, double *knodes, double k_interval){
    double kmin = knodes[0];
    double k = (double)(kmin + k_ind*k_interval);
    return k;
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  Taken from CMBfast jLgen.F and adpated to C

c  Calculates the spherical bessel function j_l(x)

c  and optionally its derivative for real x and integer l>=0.

c  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from

c  G.N.Watson, A Treatise on the Theory of Bessel Functions,
c  2nd Edition (Cambridge University Press, 1944).
c  Higher terms in expansion for x near l given by
c  Airey in Phil. Mag. 31, 520 (1916).

c  This approximation is accurate to near 0.1% at the boundaries
c  between the asymptotic regions; well away from the boundaries
c  the accuracy is better than 10^{-5}. The derivative accuracy
c  is somewhat worse than the function accuracy but still better
c  than 1%.

c  Point *jlp initially to a negative value to forego calculating
c  the derivative; point it to a positive value to do the derivative
c  also (Note: give it a definite value before the calculation
c  so it's not pointing at junk.) The derivative calculation requires

c  only arithmetic operations, plus evaluation of one sin() for the
c  x>>l region.

c  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu
c  This C version only computes j_l(x)*/

double sjl(int l, double x) {


  double jl,nu, nu2,ax,ax2,beta,beta2,beta4,beta6;
  double sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1;
  double cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b;
  double trigarg,trigcos,expterm,prefactor,llimit,ulimit,fl;

  double ROOTPI = 1.772453851;
  double GAMMA1 = 2.6789385347;           /* Gamma function of 1/3 */
  double GAMMA2 = 1.3541179394;           /* Gamma function of 2/3 */


  ax = sqrt(x*x);
  fl = (double)l;

  beta = pow(fl,0.325);
  llimit=1.31*beta;   /* limits of asymptotic regions; fitted */
  ulimit=1.48*beta;

  nu= fl + 0.5;

  nu2=nu*nu;

  if (l<0) {
    printf("Bessel function index < 0\n");
    exit(2);
  }

  /************* Use closed form for l<6 **********/

  if (l<6) {

    sx=sin(ax);
    cx=cos(ax);
    ax2=ax*ax;

    if (l==0) {
      if(ax > 0.001) {
	jl=(sx/ax);
      } else {
	jl=(1.0-ax2/6.0);
      }                       /* small x trap */
    }

    if (l==1) {
      if (ax > 0.001) {
	jl=((sx/ax -cx)/ax);
      } else {
	jl=(ax/3.0);
      }
    }

    if(l==2) {
      if(ax > 0.001) {
	jl=((-3.0*cx/ax-sx*(1.0-3.0/ax2))/ax);
      } else {
	jl=(ax2/15.0);
      }
    }

    if(l==3) {
      if(ax > 0.001) {
	jl=((cx*(1.0-15.0/ax2)-sx*(6.0-15.0/ax2)/ax)/ax);
      } else {
	jl=(ax*ax2/105.0);
      }
    }

    if(l==4) {
      if(ax > 0.001) {
        jl=((sx*(1.0-45.0/(ax*ax)+105.0/(ax*ax*ax*ax))
	     +cx*(10.0-105.0/(ax*ax))/ax)/ax);
      } else {
	jl=(ax2*ax2/945.0);
      }
    }

    if(l==5) {
      if(ax > 0.001) {
        jl=((sx*(15.0-420.0/(ax*ax)+945.0
		 /(ax*ax*ax*ax))/ax -cx*(1.0-105.0/(ax*ax)+945.0
					 /(ax*ax*ax*ax)))/ax);
      } else {
	jl=(ax2*ax2*ax/10395.0);
      }
    }

    return jl;


  }

  /********************** x=0 **********************/
  if (ax < 1.e-30) {
    jl=0.0;
    return jl;
  }

  /*************** Region 1: x << l ****************/
  if (ax <= fl+0.5-llimit) {

      //beta=acosh(nu/ax)
      if (nu/ax < 1.) printf("trouble with acosh\n");
      beta = log(nu/ax + sqrt((nu/ax)*(nu/ax) - 1.));
      //(4.6.21)
      cotb=nu/sqrt(nu*nu-ax*ax);      /* cotb=coth(beta) */
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      secb=ax/nu;
      sec2b=secb*secb;
      sec4b=sec2b*sec2b;
      sec6b=sec4b*sec2b;
      sum1=2.0+3.0*sec2b;
      expterm=sum1*cot3b/(24.0*nu);
      sum2=4.0+sec2b;
      expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2);
      sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b;
      expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2);
      sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b;
      expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2);
      expterm=exp(-nu*beta+nu/cotb-expterm);
      prefactor=sqrt(cotb/secb)/(2.0*nu);
      jl=(prefactor*expterm);

      return jl;

  }
  /**************** Region 2: x >> l ****************/
  if (ax >= fl+0.5+ulimit) {

        beta=acos(nu/ax);
        cotb=nu/sqrt(ax*ax-nu*nu);      /* cotb=cot(beta) */
        cot3b=cotb*cotb*cotb;
        cot6b=cot3b*cot3b;
        secb=ax/nu;
        sec2b=secb*secb;
        sec4b=sec2b*sec2b;
        sec6b=sec4b*sec2b;
        trigarg=nu/cotb - nu*beta - M_PI/4.0;
        sum1=2.0+3.0*sec2b;
        trigarg = trigarg - sum1*cot3b/(24.0*nu);
        sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b;
        trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2);
        trigcos=cos(trigarg);
        sum2=4.0+sec2b;
        expterm=sum2*sec2b*cot6b/(16.0*nu2);
        sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b;
        expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2);
        expterm=exp(-expterm);
        prefactor=sqrt(cotb/secb)/nu;
        jl=(prefactor*expterm*trigcos);

	return jl;

  }
  /***************** Region 3: x near l ****************/

  beta=ax-nu;

  beta2=beta*beta;
  beta4=beta2*beta2;
  beta6=beta2*beta4;
  sx=6.0/ax;
  sx2=sx*sx;
  cx=sqrt(sx);

  secb=pow(sx,0.333333333);

  sec2b=secb*secb;

  deriv1=GAMMA1*secb;
  deriv1= deriv1+ beta*GAMMA2*sec2b;
  sum1=(beta2/6.0-1.0/15.0)*beta;
  deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0;
  sum2=beta4/24.0-beta2/24.0+1.0/280.0;
  deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0;
  sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0;
  deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0;
  sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0-13.0/31500.0)*beta;
  deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0;
  sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0
	-121.0*beta2/907200.0 + 7939.0/232848000.0)*beta;
  deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0;

  jl=(deriv1*cx/(12.0*ROOTPI));

  return jl;
}
