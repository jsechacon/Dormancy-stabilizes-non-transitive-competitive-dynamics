
###################################################################################
########### EULER-MURAYAMA SCHEME FOR SIMULATING THE SCALING LIMIT ###############
###################################################################################

#include <Rcpp.h>
using namespace Rcpp;

// Helper: project a vector to the simplex
inline void project_simplex(NumericVector &v) {
    int n = v.size();
    double s = 0.0;

    for (int i = 0; i < n; i++) {
        if (v[i] < 0) v[i] = 0;
        s += v[i];
    }
    if (s > 1e-15) {
        for (int i = 0; i < n; i++) v[i] /= s;
    }
}

// [[Rcpp::export]]
List simulate_wf_seedbank_cpp(
        double T,
        double dt,
        int noTypes,
        NumericVector x0,
        NumericVector y0,
        NumericVector drift_const,   // dummy (ignored, see R wrapper)
        Function muR,                // drift from R
        Function SigmaR,             // diffusion matrix from R
        double c1,
        double c2,
        double sigma
) {
    int nSteps = (int) std::ceil(T / dt);

    NumericVector x = clone(x0);
    NumericVector y = clone(y0);
    NumericVector dx(noTypes);
    NumericVector dy(noTypes);
    NumericVector dB(noTypes);

    double time = 0.0;

    for (int step = 0; step < nSteps; step++) {

        // --- drift ---
        dx = as<NumericVector>(muR(x));
        dy = c2 * (x - y);

        for (int i = 0; i < noTypes; i++) {
            dx[i] = dx[i] * dt + c1 * (y[i] - x[i]) * dt;
            dy[i] = dy[i] * dt;
        }

        // --- Brownian increment ---
        for (int i = 0; i < noTypes; i++) {
            dB[i] = R::rnorm(0.0, std::sqrt(dt));
        }

        // --- diffusion ---
        NumericMatrix S = as<NumericMatrix>(SigmaR(x));
        NumericVector diff(noTypes);
        for (int i = 0; i < noTypes; i++) {
            double sum = 0.0;
            for (int j = 0; j < noTypes; j++) {
                sum += S(i, j) * dB[j];
            }
            diff[i] = std::sqrt(sigma) * sum;
        }

        // --- apply updates ---
        for (int i = 0; i < noTypes; i++) {
            x[i] += dx[i] + diff[i];
            y[i] += dy[i];
        }

        // --- project to simplexes ---
        project_simplex(x);
        project_simplex(y);

        time += dt;

        // fixation check
        double xmax = 0.0;
        int idx = 0;
        for (int i = 0; i < noTypes; i++) {
            if (x[i] > xmax) {
                xmax = x[i];
                idx = i + 1;
            }
        }
        if (xmax == 1.0) {
            return List::create(
                _["FixTime"] = time,
                _["FixType"] = idx
            );
        }
    }

    // failed to fix
    return List::create(
        _["FixTime"] = -1.0,
        _["FixType"] = -1
    );
}
