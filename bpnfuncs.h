#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

void binomial(mpz_t out,int n,int k);
void inpdf(mpfr_t *out,mpz_t n,int psize,double pini,double pinc,unsigned long mu,mpfr_t bcs[],mpfr_prec_t prec);
void hypergeometricpdf(mpfr_t *out,mpz_t mu, mpz_t d, unsigned long int mi, unsigned long int dini, unsigned long int dmax, mpz_t bc_b_ai[], mpz_t mu_bc[]);
void pphi(mpfr_t *pout,mpz_t mu,mpz_t dmp,int m,int dini, int dmax,mpz_t bc_b_ai[], mpz_t mu_bc[], mpfr_prec_t prec);
void uniquecodes(mpfr_t *ucodes, int psiz, mpz_t n, unsigned long int mu, mpz_t ncodes ,mpfr_t *pdf,mpfr_prec_t prec);
void scjoint(mpfr_t *distptr, int mu, int gamma, int psiz, int dini, int dmax, mpfr_t *bcs,mpfr_t *binpdf,mpfr_t *pp,mpfr_prec_t prec);
void tabulate_bins_fr(mpfr_t *bcs, int Nini, int Nend, int mini, int mend, mpfr_prec_t prec);
void tabulate_bins_z(mpz_t bcs[], int Nini, int Nend, int mini, int mend);
void net_entropy(mpfr_t *ent, mpfr_t *ipdf, mpfr_t *jdist, mpfr_t *pp, mpfr_t *gamma_bcs, mpfr_t *p0, int psiz, int dini, int dmax, int mu, int gamma, mpfr_prec_t prec);
void patt_sep_ou(mpfr_t *psep, mpfr_t *ipdf, mpfr_t *pp, mpfr_t *jpdf, mpfr_t *gamma_bcs, int psiz, int dini, int dmax, int mu, int gamma, mpfr_prec_t prec);
void patt_sep_in(mpfr_t *psep, mpfr_t *ipdf, mpfr_t *mu_bcs, int psiz, int mu, double pinc, mpfr_prec_t prec);

