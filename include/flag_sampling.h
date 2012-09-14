// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_SAMPLING
#define FLAG_SAMPLING

/* FOR FUTURE IMPROVEMENTS // multi-scheme support
enum ssht_methods { 
	MW, 
	MWSS, 
	GL, 
	DH 
};
*/

/*!
 * Allocate FLAG sampling.
 *
 * \param[out]  rs Radial coordinates.
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[out]  laguweights Laguerre radial weights for FLAG transform.
 * \param[in]  R Radial boundary / limit.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N);

/*!
 * Allocate SSHT sampling.
 *
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void ssht_allocate_sampling(double **thetas, double **phis, int L);

/*!
 * Compute SSHT MW sampling.
 *
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void ssht_sampling_mw(double *thetas, double *phis, int L);

/*!
 * Compute FLAG sampling.
 *
 * \param[out]  rs Radial coordinates.
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[out]  laguweights Laguerre radial weights for FLAG transform.
 * \param[in]  R Radial boundary / limit.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N);

#endif