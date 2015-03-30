

#include "fitting.h"
extern "C" {
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
}
#include <stdint.h>


void fitPlane(const std::vector<Vector> & x, float & a, float & b, float & d)
{
    gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(x.size(), 3);
    gsl_matrix *X = gsl_matrix_calloc(x.size(), 3);
    gsl_vector *Y = gsl_vector_alloc(x.size());
    gsl_vector *beta = gsl_vector_alloc(3);
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);
    double chisq;
    for(uint32_t i = 0;i < x.size();i++){
        gsl_vector_set(Y, i, x[i].z());
        gsl_matrix_set(X, i, 0, 1);
        gsl_matrix_set(X, i, 1, x[i].x());
        gsl_matrix_set(X, i, 2, x[i].y());
    }
    gsl_multifit_linear(X, Y, beta, cov, &chisq, ws);
    /* cz = d + ax + by */
    d = gsl_vector_get(beta, 0);
    a = gsl_vector_get(beta, 1);
    b = gsl_vector_get(beta, 2);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(Y);
    gsl_vector_free(beta);
    gsl_multifit_linear_free(ws);
}

