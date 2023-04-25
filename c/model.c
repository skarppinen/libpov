#include <math.h>
#include <assert.h>
#include "random-sweep.h"

// Get index of (i, j, k) from array A with slices of size `m` * `n`,
// that is, m rows and n columns.
// Slices are the last dimension.
// The naming comes from Column Major Element: Row, Column, Slice
#define CMERCS(i, j, k, m, n) ((k) * (m) * (n) + (j) * (m) + (i)) 

typedef struct {
    unsigned int nS; // Number of stands.
    unsigned int nT; // Number of time points.
    unsigned int nA; // Number of assortments.
    unsigned int nI; // Number of different inventory methods.
    double * muprior; // Prior volumes, `nS` x `nA` array.
    double * sigma2prior; // Prior volume variances, `nS` x `nA` array.
    double * sigma2meas; // Volume measurement variances, `nI` x `nS` x `nA` array.
    double * demands; // Demands per time and assortment, `nT` x `nA` array.
} StandData;

typedef struct {
    unsigned int nS; // Number of stands.
    unsigned int nA; // Number of assortments.
    double * muplus; // Posterior mean volumes, `nS` x `nA` matrix.
    double * sigma2plus; // Posterior volume variances, `nS` x `nA` matrix.
    double * y; // Temporary for simulated observations, `nS` x `nA` matrix.
} VolumePosterior;

typedef struct {
    unsigned long * arr; // A vector of length `n` with random seeds.
    unsigned int n; // Number of random seeds.
} SeedVector;

void print_vp(const VolumePosterior * vp) {
    printf("VOLUME POSTERIOR\n");
    printf("nS = %u, nA = %u\n", vp->nS, vp->nA);
    printf("First three values:\n");
    printf("sigma2plus: %lf, %lf, %lf\n", vp->sigma2plus[0], vp->sigma2plus[1], vp->sigma2plus[2]);
    printf("muplus: %lf, %lf, %lf\n", vp->muplus[0], vp->muplus[1], vp->muplus[2]);
    printf("y: %lf, %lf, %lf\n", vp->y[0], vp->y[1], vp->y[2]);
}

static void set_to_zero(double * x, const unsigned int n) {
    for (unsigned int i = 0; i < n; i++)
        x[i] = 0.0;
}

static void u_to_l_tri(double * cm_mat, const unsigned int dim) {
	unsigned int step = 0;
	unsigned int colstart = 0;
	for (unsigned int j = 0; j < dim - 1; ++j) {
		step = dim - 1;
		colstart = j * dim;
		for (unsigned int i = j + 1; i < dim; ++i) {
			cm_mat[colstart + i] = cm_mat[colstart + i + step];
			step += dim - 1;
		}
	}
}

void build_problem(RandomSweepStorage *restrict rs, 
                   const VolumePosterior *restrict vp,
                   const double *restrict demands) {
    double r = 0.0;
    for (unsigned int i = 0; i < vp->nA * rs->nT; i++) {
        r += demands[i] * demands[i]; 
    }
    rs->r = r;

    // Computes `rs.C = -2.0 * muplus * transpose(demands)`.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                vp->nS, rs->nT, vp->nA, 
                -2.0, vp->muplus, vp->nS, 
                demands, rs->nT,
                0.0, rs->C, vp->nS);

    // Computes `rs.Qb`. 
    set_to_zero(rs->Qb, vp->nS * vp->nS);
    for (unsigned int a = 0; a < vp->nA; a++) {
        // Here using striding by `nS + 1` to update diagonal elements of Qb.
        // This adds 2 * Sigma^a to Qb. (Sigma^a is diagonal)
        cblas_daxpy(vp->nS, 
                    2.0, vp->sigma2plus + a * vp->nS, 1,
                    rs->Qb, vp->nS + 1); 
    }
    // This adds sum_{a} 2 * muplus^a * transpose(muplus^a) to Qb.
    // syrk only writes to the upper triangle of Qb here.
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                vp->nS, vp->nA, 
                2.0, vp->muplus, vp->nS,
                1.0, rs->Qb, vp->nS);
    u_to_l_tri(rs->Qb, vp->nS); // Fill lower triangle with what is on upper triangle.
};

static void normal_posterior_mean_var(double * mupost, double * sigma2post,
                                      double y, double muprior, double sigma2prior,
                                      double sigma2meas) { 
    *sigma2post = 1.0 / ((1.0 / sigma2prior) + (1.0 / sigma2meas));
    *mupost = *sigma2post * (muprior / sigma2prior + y / sigma2meas);
}

void simulate_data(double *restrict y, const unsigned int *restrict xI, 
                   const StandData *restrict stda, gsl_rng *restrict rng) {
    double priorm = 0.0;
    double priorv = 0.0;
    double measv = 0.0;
    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            measv = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            priorm = stda->muprior[CMERC(s, a, stda->nS)];
            priorv = stda->sigma2prior[CMERC(s, a, stda->nS)];
            y[CMERC(s, a, stda->nS)] = priorm + gsl_ran_gaussian_ziggurat(rng, sqrt(priorv + measv)); 
        }
    }
}

void compute_posterior(VolumePosterior *restrict vp,
                       const unsigned int *restrict xI,
                       const StandData *restrict stda) {
    double priorm = 0.0;
    double priorv = 0.0;
    double measv = 0.0;
    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            measv = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            priorm = stda->muprior[CMERC(s, a, stda->nS)];
            priorv = stda->sigma2prior[CMERC(s, a, stda->nS)];
            normal_posterior_mean_var(vp->muplus + CMERC(s, a, vp->nS), 
                                      vp->sigma2plus + CMERC(s, a, vp->nS),
                                      vp->y[CMERC(s, a, vp->nS)], priorm, priorv, measv);
        }
    }
}



void random_posterior(VolumePosterior *restrict vp,
                      const unsigned int *restrict xI,
                      const StandData *restrict stda,
                      gsl_rng *restrict rng) {
    simulate_data(vp->y, xI, stda, rng);
    compute_posterior(vp, xI, stda);
}

double PoV(const unsigned int *restrict xI,
           RandomSweepStorage *restrict rs,
           VolumePosterior *restrict vp,
           const StandData *restrict stda, 
           const SeedVector *restrict seeds,  
           const unsigned int maxsweeps,
           const unsigned int ninits) {

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus); 
    double out = 0.0; 
    double result = 0.0;
    for (unsigned int n = 1; n <= seeds->n; n++) {
        gsl_rng_set(rng, seeds->arr[n - 1]);
        random_posterior(vp, xI, stda, rng);
        build_problem(rs, vp, stda->demands); 
        random_sweep(&result, rs, maxsweeps, ninits, rng); 
        out = out * (n - 1) / n + result / n; // Running mean. 
    }
    gsl_rng_free(rng);
    return out;
}


