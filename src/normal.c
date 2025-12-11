/* normal.c  - core C implementation for Normal distribution
 *
 * To compile with Rtools on Windows:
 *   R CMD SHLIB normal.c
 *
 * This will create normal.dll, which can be loaded from R and
 * called via .C("normal_pdf_c", ...) etc.
 */

#include <math.h>

/* constant */
#ifndef PI
#define PI 3.14159265358979323846
#endif

/* -------------------------
   Internal scalar helpers
   ------------------------- */

/* Normal PDF */
static double normal_pdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * PI));
}

/* Normal CDF using an approximation to the error function
   Approximation from Abramowitz & Stegun, 7.1.26
*/
static double normal_cdf(double x, double mu, double sigma) {
    double z = (x - mu) / (sigma * sqrt(2.0));
    double t, erf_approx, sign;

    sign = (z < 0.0) ? -1.0 : 1.0;
    z = fabs(z);

    /* erf approximation parameters */
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p  = 0.3275911;

    t = 1.0 / (1.0 + p * z);
    erf_approx = 1.0 - (((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t) * exp(-z*z);

    return 0.5 * (1.0 + sign * erf_approx);
}

/* -------------------------
   R .C interface (vectorised)
   All arguments are pointers
   ------------------------- */

/* pdf: compute pdf for vector x of length n
   results written into out[i]
*/
void normal_pdf_c(int *n, double *x,
                  double *mu, double *sigma,
                  double *out) {
    int i;
    for (i = 0; i < *n; i++) {
        out[i] = normal_pdf(x[i], *mu, *sigma);
    }
}

/* cdf: compute cdf for vector x of length n */
void normal_cdf_c(int *n, double *x,
                  double *mu, double *sigma,
                  double *out) {
    int i;
    for (i = 0; i < *n; i++) {
        out[i] = normal_cdf(x[i], *mu, *sigma);
    }
}

/* prob: compute interval probability P(lower[i] <= X <= upper[i])
   lower / upper are vectors of length n
*/
void normal_prob_c(int *n,
                   double *lower, double *upper,
                   double *mu, double *sigma,
                   double *out) {
    int i;
    for (i = 0; i < *n; i++) {
        double cdf_upper = normal_cdf(upper[i], *mu, *sigma);
        double cdf_lower = normal_cdf(lower[i], *mu, *sigma);
        out[i] = cdf_upper - cdf_lower;
    }
}
