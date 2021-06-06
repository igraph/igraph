#include <igraph.h>

#include "bench.h"
#include "../../vendor/plfit/sampling.h"

igraph_vector_t data;

inline double rpareto(double xmin, double alpha) {
    /* 1-u is used in the base here because we want to avoid the case of
     * sampling zero */
    return pow(1 - RNG_UNIF01(), -1.0 / alpha) * xmin;
}

inline double rzeta(long int xmin, double alpha) {
    double u, v, t;
    long int x;
    double alpha_minus_1 = alpha-1;
    double minus_1_over_alpha_minus_1 = -1.0 / (alpha-1);
    double b;
    double one_over_b_minus_1;

    xmin = (long int) round(xmin);

    /* Rejection sampling for the win. We use Y=floor(U^{-1/alpha} * xmin) as the
     * envelope distribution, similarly to Chapter X.6 of Luc Devroye's book
     * (where xmin is assumed to be 1): http://luc.devroye.org/chapter_ten.pdf
     *
     * Some notes that should help me recover what I was doing:
     *
     * p_i = 1/zeta(alpha, xmin) * i^-alpha
     * q_i = (xmin/i)^{alpha-1} - (xmin/(i+1))^{alpha-1}
     *     = (i/xmin)^{1-alpha} - ((i+1)/xmin)^{1-alpha}
     *     = [i^{1-alpha} - (i+1)^{1-alpha}] / xmin^{1-alpha}
     *
     * p_i / q_i attains its maximum at xmin=i, so the rejection constant is:
     *
     * c = p_xmin / q_xmin
     *
     * We have to accept the sample if V <= (p_i / q_i) * (q_xmin / p_xmin) =
     * (i/xmin)^-alpha * [xmin^{1-alpha} - (xmin+1)^{1-alpha}] / [i^{1-alpha} - (i+1)^{1-alpha}] =
     * [xmin - xmin^alpha / (xmin+1)^{alpha-1}] / [i - i^alpha / (i+1)^{alpha-1}] =
     * xmin/i * [1-(xmin/(xmin+1))^{alpha-1}]/[1-(i/(i+1))^{alpha-1}]
     *
     * In other words (and substituting i with X, which is the same),
     *
     * V * (X/xmin) <= [1 - (1+1/xmin)^{1-alpha}] / [1 - (1+1/i)^{1-alpha}]
     *
     * Let b := (1+1/xmin)^{alpha-1} and let T := (1+1/i)^{alpha-1}. Then:
     *
     * V * (X/xmin) <= [(b-1)/b] / [(T-1)/T]
     * V * (X/xmin) * (T-1) / (b-1) <= T / b
     *
     * which is the same as in Devroye's book, except for the X/xmin term, and
     * the definition of b.
     */
    b = pow(1 + 1.0/xmin, alpha_minus_1);
    one_over_b_minus_1 = 1.0/(b-1);
    do {
        do {
            u = RNG_UNIF01();
            v = RNG_UNIF01();
            /* 1-u is used in the base here because we want to avoid the case of
             * having zero in x */
            x = (long int) floor(pow(1-u, minus_1_over_alpha_minus_1) * xmin);
        } while (x < xmin);
        t = pow((x+1.0)/x, alpha_minus_1);
    } while (v*x*(t-1)*one_over_b_minus_1*b > t*xmin);

    return x;
}

int generate_continuous(double xmin, double alpha, size_t num_samples) {
    IGRAPH_CHECK(igraph_vector_resize(&data, num_samples));

    RNG_BEGIN();
    for (size_t i = 0; i < num_samples; i++) {
        VECTOR(data)[i] = rpareto(xmin, alpha);
    }
    RNG_END();

    return IGRAPH_SUCCESS;
}

int generate_discrete(double xmin, double alpha, size_t num_samples) {
    IGRAPH_CHECK(igraph_vector_resize(&data, num_samples));

    RNG_BEGIN();
    for (size_t i = 0; i < num_samples; i++) {
        VECTOR(data)[i] = rzeta(xmin, alpha);
    }
    RNG_END();

    return IGRAPH_SUCCESS;
}

int fit_continuous(double known_xmin) {
    igraph_plfit_result_t result;
    IGRAPH_CHECK(igraph_power_law_fit(&data, &result, known_xmin, /* force_continuous = */ 1));
    return IGRAPH_SUCCESS;
}

int fit_discrete(double known_xmin) {
    igraph_plfit_result_t result;
    IGRAPH_CHECK(igraph_power_law_fit(&data, &result, known_xmin, 0));
    return IGRAPH_SUCCESS;
}

int main() {
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&data, 0);

    BENCH_INIT();

    generate_continuous(1, 3, 100000);
    BENCH(" 1 Continuous, xmin = 1, alpha = 3, samples = 100K, fitting alpha only", fit_continuous(1));

    generate_continuous(1, 3, 200000);
    BENCH(" 2 Continuous, xmin = 1, alpha = 3, samples = 200K, fitting alpha only", fit_continuous(1));

    generate_continuous(1, 3, 500000);
    BENCH(" 3 Continuous, xmin = 1, alpha = 3, samples = 500K, fitting alpha only", fit_continuous(1));

    generate_continuous(1, 3, 5000);
    BENCH(" 4 Continuous, xmin = 1, alpha = 3, samples = 5K, fitting xmin and alpha", fit_continuous(-1));

    generate_continuous(1, 3, 10000);
    BENCH(" 5 Continuous, xmin = 1, alpha = 3, samples = 10K, fitting xmin and alpha", fit_continuous(-1));

    generate_continuous(1, 3, 15000);
    BENCH(" 6 Continuous, xmin = 1, alpha = 3, samples = 15K, fitting xmin and alpha", fit_continuous(-1));

    generate_discrete(3, 3, 1000000);
    BENCH(" 7 Discrete, xmin = 3, alpha = 3, samples = 1M, fitting alpha only", fit_discrete(3));

    generate_discrete(3, 3, 5000000);
    BENCH(" 8 Discrete, xmin = 3, alpha = 3, samples = 5M, fitting alpha only", fit_discrete(3));

    generate_discrete(3, 3, 10000000);
    BENCH(" 9 Discrete, xmin = 3, alpha = 3, samples = 10M, fitting alpha only", fit_discrete(3));

    generate_discrete(3, 3, 1000000);
    BENCH("10 Discrete, xmin = 3, alpha = 3, samples = 1M, fitting xmin and alpha", fit_discrete(-1));

    generate_discrete(3, 3, 5000000);
    BENCH("11 Discrete, xmin = 3, alpha = 3, samples = 5M, fitting xmin and alpha", fit_discrete(-1));

    generate_discrete(3, 3, 10000000);
    BENCH("12 Discrete, xmin = 3, alpha = 3, samples = 10M, fitting xmin and alpha", fit_discrete(-1));

    igraph_vector_destroy(&data);

    return 0;
}
