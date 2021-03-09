#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <stdexcept>
using namespace prpack;
using namespace std;

void prpack_solver::initialize() {
    geg = NULL;
    gsg = NULL;
    sg = NULL;
    sccg = NULL;
	owns_bg = true;
}

prpack_solver::prpack_solver(const prpack_csc* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(const prpack_int64_csc* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(const prpack_csr* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(const prpack_edge_list* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(prpack_base_graph* g, bool owns_bg) {
    initialize();
	this->owns_bg = owns_bg;
    TIME(read_time, bg = g);
}

prpack_solver::prpack_solver(const char* filename, const char* format, const bool weighted) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(filename, format, weighted));
}

prpack_solver::~prpack_solver() {
	if (owns_bg) {
		delete bg;
	}
    delete geg;
    delete gsg;
    delete sg;
    delete sccg;
}

int prpack_solver::get_num_vs() {
    return bg->num_vs;
}

prpack_result* prpack_solver::solve(const double alpha, const double tol, const char* method) {
    return solve(alpha, tol, NULL, NULL, method);
}

prpack_result* prpack_solver::solve(
        const double alpha,
        const double tol,
        const double* u,
        const double* v,
        const char* method) {
    double preprocess_time = 0;
    double compute_time = 0;
    prpack_result* ret = NULL;
    // decide which method to run
    string m;
    if (strcmp(method, "") != 0)
        m = string(method);
    else {
        if (bg->num_vs < 128)
            m = "ge";
        else if (sccg != NULL)
            m = "sccgs";
        else if (sg != NULL)
            m = "sg";
        else
            m = "sccgs";
        if (u != v)
            m += "_uv";
    }
    // run the appropriate method
    if (m == "ge") {
        if (geg == NULL) {
            TIME(preprocess_time, geg = new prpack_preprocessed_ge_graph(bg));
        }
        TIME(compute_time, ret = solve_via_ge(
                alpha,
                tol,
                geg->num_vs,
                geg->matrix,
                u));
    } else if (m == "ge_uv") {
        if (geg == NULL) {
            TIME(preprocess_time, geg = new prpack_preprocessed_ge_graph(bg));
        }
        TIME(compute_time, ret = solve_via_ge_uv(
                alpha,
                tol,
                geg->num_vs,
                geg->matrix,
                geg->d,
                u,
                v));
    } else if (m == "gs") {
        if (gsg == NULL) {
            TIME(preprocess_time, gsg = new prpack_preprocessed_gs_graph(bg));
        }
        TIME(compute_time, ret = solve_via_gs(
                alpha,
                tol,
                gsg->num_vs,
                gsg->num_es,
                gsg->heads,
                gsg->tails,
                gsg->vals,
                gsg->ii,
                gsg->d,
                gsg->num_outlinks,
                u,
                v));
    } else if (m == "gserr") {
        if (gsg == NULL) {
            TIME(preprocess_time, gsg = new prpack_preprocessed_gs_graph(bg));
        }
        TIME(compute_time, ret = solve_via_gs_err(
                alpha,
                tol,
                gsg->num_vs,
                gsg->num_es,
                gsg->heads,
                gsg->tails,
                gsg->ii,
                gsg->num_outlinks,
                u,
                v));
    } else if (m == "sgs") {
        if (sg == NULL) {
            TIME(preprocess_time, sg = new prpack_preprocessed_schur_graph(bg));
        }
        TIME(compute_time, ret = solve_via_schur_gs(
                alpha,
                tol,
                sg->num_vs,
                sg->num_no_in_vs,
                sg->num_no_out_vs,
                sg->num_es,
                sg->heads,
                sg->tails,
                sg->vals,
                sg->ii,
                sg->d,
                sg->num_outlinks,
                u,
                sg->encoding,
                sg->decoding));
    } else if (m == "sgs_uv") {
        if (sg == NULL) {
            TIME(preprocess_time, sg = new prpack_preprocessed_schur_graph(bg));
        }
        TIME(compute_time, ret = solve_via_schur_gs_uv(
                alpha,
                tol,
                sg->num_vs,
                sg->num_no_in_vs,
                sg->num_no_out_vs,
                sg->num_es,
                sg->heads,
                sg->tails,
                sg->vals,
                sg->ii,
                sg->d,
                sg->num_outlinks,
                u,
                v,
                sg->encoding,
                sg->decoding));
    } else if (m == "sccgs") {
        if (sccg == NULL) {
            TIME(preprocess_time, sccg = new prpack_preprocessed_scc_graph(bg));
        }
        TIME(compute_time, ret = solve_via_scc_gs(
                alpha,
                tol,
                sccg->num_vs,
                sccg->num_es_inside,
                sccg->heads_inside,
                sccg->tails_inside,
                sccg->vals_inside,
                sccg->num_es_outside,
                sccg->heads_outside,
                sccg->tails_outside,
                sccg->vals_outside,
                sccg->ii,
                sccg->d,
                sccg->num_outlinks,
                u,
                sccg->num_comps,
                sccg->divisions,
                sccg->encoding,
                sccg->decoding));
    } else if (m == "sccgs_uv") {
        if (sccg == NULL) {
            TIME(preprocess_time, sccg = new prpack_preprocessed_scc_graph(bg));
        }
        TIME(compute_time, ret = solve_via_scc_gs_uv(
                alpha,
                tol,
                sccg->num_vs,
                sccg->num_es_inside,
                sccg->heads_inside,
                sccg->tails_inside,
                sccg->vals_inside,
                sccg->num_es_outside,
                sccg->heads_outside,
                sccg->tails_outside,
                sccg->vals_outside,
                sccg->ii,
                sccg->d,
                sccg->num_outlinks,
                u,
                v,
                sccg->num_comps,
                sccg->divisions,
                sccg->encoding,
                sccg->decoding));
    } else {
        throw invalid_argument("Unknown method specified for PRPACK: '" + m + "'.");
    }
    ret->method = m;
    ret->read_time = read_time;
    ret->preprocess_time = preprocess_time;
    ret->compute_time = compute_time;
    ret->num_vs = bg->num_vs;
    ret->num_es = bg->num_es;
    return ret;
}

// VARIOUS SOLVING METHODS ////////////////////////////////////////////////////////////////////////

prpack_result* prpack_solver::solve_via_ge(
        const double alpha,
        const double tol,
        const int num_vs,
        const double* matrix,
        const double* uv) {
    prpack_result* ret = new prpack_result();
    // initialize uv values
    const double uv_const = 1.0/num_vs;
    const int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? uv : &uv_const;
    // create matrix A
    double* A = new double[num_vs*num_vs];
    for (int i = 0; i < num_vs*num_vs; ++i)
        A[i] = -alpha*matrix[i];
    for (int i = 0; i < num_vs*num_vs; i += num_vs + 1)
        ++A[i];
    // create vector b
    double* b = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        b[i] = uv[uv_exists*i];
    // solve and normalize
    ge(num_vs, A, b);
    normalize(num_vs, b);
    // clean up and return
    delete[] A;
    ret->num_es_touched = -1;
    ret->x = b;
    return ret;
}

prpack_result* prpack_solver::solve_via_ge_uv(
        const double alpha,
        const double tol,
        const int num_vs,
        const double* matrix,
        const double* d,
        const double* u,
        const double* v) {
    prpack_result* ret = new prpack_result();
    // initialize u and v values
    const double u_const = 1.0/num_vs;
    const double v_const = 1.0/num_vs;
    const int u_exists = (u) ? 1 : 0;
    const int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // create matrix A
    double* A = new double[num_vs*num_vs];
    for (int i = 0; i < num_vs*num_vs; ++i)
        A[i] = -alpha*matrix[i];
    for (int i = 0, inum_vs = 0; i < num_vs; ++i, inum_vs += num_vs)
        for (int j = 0; j < num_vs; ++j)
            A[inum_vs + j] -= alpha*u[u_exists*i]*d[j];
    for (int i = 0; i < num_vs*num_vs; i += num_vs + 1)
        ++A[i];
    // create vector b
    double* b = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        b[i] = (1 - alpha)*v[v_exists*i];
    // solve
    ge(num_vs, A, b);
    // clean up and return
    delete[] A;
    ret->num_es_touched = -1;
    ret->x = b;
    return ret;
}

// Vanilla Gauss-Seidel.
prpack_result* prpack_solver::solve_via_gs(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es,
        const int* heads,
        const int* tails,
        const double* vals,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* u,
        const double* v) {
    prpack_result* ret = new prpack_result();
    const bool weighted = vals != NULL;
    // initialize u and v values
    const double u_const = 1.0/num_vs;
    const double v_const = 1.0/num_vs;
    const int u_exists = (u) ? 1 : 0;
    const int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        x[i] = 0;
    // initialize delta
    double delta = 0;
    // run Gauss-Seidel
    ret->num_es_touched = 0;
    double err = 1, c = 0;
    do {
        if (weighted) {
            for (int i = 0; i < num_vs; ++i) {
                double new_val = 0;
                const int start_j = tails[i];
                const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
                for (int j = start_j; j < end_j; ++j)
                    // TODO: might want to use compensation summation for large: end_j - start_j
                    new_val += x[heads[j]]*vals[j];
                new_val = alpha*new_val + (1 - alpha)*v[v_exists*i];
                delta -= alpha*x[i]*d[i];
                new_val += delta*u[u_exists*i];
                new_val /= 1 - alpha*(d[i]*u[u_exists*i] + (1 - d[i])*ii[i]);
                delta += alpha*new_val*d[i];
                COMPENSATED_SUM(err, x[i] - new_val, c);
                x[i] = new_val;
            }
        } else {
            for (int i = 0; i < num_vs; ++i) {
                const double old_val = x[i]*num_outlinks[i];
                double new_val = 0;
                const int start_j = tails[i];
                const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
                for (int j = start_j; j < end_j; ++j)
                    // TODO: might want to use compensation summation for large: end_j - start_j
                    new_val += x[heads[j]];
                new_val = alpha*new_val + (1 - alpha)*v[v_exists*i];
                if (num_outlinks[i] < 0) {
                    delta -= alpha*old_val;
                    new_val += delta*u[u_exists*i];
                    new_val /= 1 - alpha*u[u_exists*i];
                    delta += alpha*new_val;
                } else {
                    new_val += delta*u[u_exists*i];
                    new_val /= 1 - alpha*ii[i];
                }
                COMPENSATED_SUM(err, old_val - new_val, c);
                x[i] = new_val/num_outlinks[i];
            }
        }
        // update iteration index
        ret->num_es_touched += num_es;
    } while (err >= tol);
    // undo num_outlinks transformation
    if (!weighted)
        for (int i = 0; i < num_vs; ++i)
            x[i] *= num_outlinks[i];
    // return results
    ret->x = x;
    return ret;
}

// Implement a gauss-seidel-like process with a strict error bound
// we return a solution with 1-norm error less than tol.
prpack_result* prpack_solver::solve_via_gs_err(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es,
        const int* heads,
        const int* tails,
        const double* ii,
        const double* num_outlinks,
        const double* u,
        const double* v) {
    prpack_result* ret = new prpack_result();
    // initialize u and v values
    const double u_const = 1.0/num_vs;
    const double v_const = 1.0/num_vs;
    const int u_exists = (u) ? 1 : 0;
    const int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // Note to Dave, we can't rescale v because we could be running this
    // same routine from multiple threads.
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i) {
        x[i] = 0.;
    }
    // initialize delta
    double delta = 0.;
    // run Gauss-Seidel, note that we store x/deg[i] throughout this 
    // iteration.
    int64_t maxedges = (int64_t)((double)num_es*std::min(
                            log(tol)/log(alpha),
                            (double)PRPACK_SOLVER_MAX_ITERS));
    ret->num_es_touched = 0;
    double err=1., c = 0.;
    do {
        // iterate through vertices
        for (int i = 0; i < num_vs; ++i) {
            double old_val = x[i]*num_outlinks[i]; // adjust back to the "true" value.
            double new_val = 0.;
            int start_j = tails[i], end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
            for (int j = start_j; j < end_j; ++j) {
                // TODO: might want to use compensation summation for large: end_j - start_j
                new_val += x[heads[j]];
            }
            new_val = alpha*new_val + alpha*ii[i]*old_val + (1.0-alpha)*v[v_exists*i];
            new_val += delta*u[u_exists*i]; // add the dangling node adjustment
            if (num_outlinks[i] < 0) {
                delta += alpha*(new_val - old_val);
            } 
            // note that new_val > old_val, but the fabs is just for 
            COMPENSATED_SUM(err, -(new_val - old_val), c);
            x[i] = new_val/num_outlinks[i];
        }
        // update iteration index
        ret->num_es_touched += num_es;
    } while (err >= tol && ret->num_es_touched < maxedges);
    if (err >= tol) {
        ret->converged = 0;
    } else {
        ret->converged = 1;
    }
    // undo num_outlinks transformation
    for (int i = 0; i < num_vs; ++i)
        x[i] *= num_outlinks[i];
    // return results
    ret->x = x;
    return ret;
}

// Gauss-Seidel using the Schur complement to separate dangling nodes.
prpack_result* prpack_solver::solve_via_schur_gs(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_no_in_vs,
        const int num_no_out_vs,
        const int num_es,
        const int* heads,
        const int* tails,
        const double* vals,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* uv,
        const int* encoding,
        const int* decoding,
        const bool should_normalize) {
    prpack_result* ret = new prpack_result();
    const bool weighted = vals != NULL;
    // initialize uv values
    const double uv_const = 1.0/num_vs;
    const int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? prpack_utils::permute(num_vs, uv, encoding) : &uv_const;
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs - num_no_out_vs; ++i)
        x[i] = uv[uv_exists*i]/(1 - alpha*ii[i])/((weighted) ? 1 : num_outlinks[i]);
    // run Gauss-Seidel for the top left part of (I - alpha*P)*x = uv
    ret->num_es_touched = 0;
    double err, c;
    do {
        // iterate through vertices
        int num_es_touched = 0;
        err = c = 0;
        #pragma omp parallel for firstprivate(c) reduction(+:err, num_es_touched) schedule(dynamic, 64)
        for (int i = num_no_in_vs; i < num_vs - num_no_out_vs; ++i) {
            double new_val = 0;
            const int start_j = tails[i];
            const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
            if (weighted) {
                for (int j = start_j; j < end_j; ++j)
                    // TODO: might want to use compensation summation for large: end_j - start_j
                    new_val += x[heads[j]]*vals[j];
                COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]), c);
                new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                x[i] = new_val;
            } else {
                for (int j = start_j; j < end_j; ++j)
                    // TODO: might want to use compensation summation for large: end_j - start_j
                    new_val += x[heads[j]];
                COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]*num_outlinks[i]), c);
                new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                x[i] = new_val/num_outlinks[i];
            }
            num_es_touched += end_j - start_j;
        }
        // update iteration index
        ret->num_es_touched += num_es_touched;
    } while (err/(1 - alpha) >= tol);
    // solve for the dangling nodes
    int num_es_touched = 0;
    #pragma omp parallel for reduction(+:num_es_touched) schedule(dynamic, 64)
    for (int i = num_vs - num_no_out_vs; i < num_vs; ++i) {
        x[i] = 0;
        const int start_j = tails[i];
        const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
        for (int j = start_j; j < end_j; ++j)
            x[i] += x[heads[j]]*((weighted) ? vals[j] : 1);
        x[i] = (alpha*x[i] + uv[uv_exists*i])/(1 - alpha*ii[i]);
        num_es_touched += end_j - start_j;
    }
    ret->num_es_touched += num_es_touched;
    // undo num_outlinks transformation
    if (!weighted)
        for (int i = 0; i < num_vs - num_no_out_vs; ++i)
            x[i] *= num_outlinks[i];
    // normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
    if (should_normalize)
        normalize(num_vs, x);
    // return results
    ret->x = prpack_utils::permute(num_vs, x, decoding);
    delete[] x;
    if (uv_exists)
        delete[] uv;
    return ret;
}

prpack_result* prpack_solver::solve_via_schur_gs_uv(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_no_in_vs,
        const int num_no_out_vs,
        const int num_es,
        const int* heads,
        const int* tails,
        const double* vals,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* u,
        const double* v,
        const int* encoding,
        const int* decoding) {
    // solve uv = u
    prpack_result* ret_u = solve_via_schur_gs(
            alpha,
            tol,
            num_vs,
            num_no_in_vs,
            num_no_out_vs,
            num_es,
            heads,
            tails,
            vals,
            ii,
            d,
            num_outlinks,
            u,
            encoding,
            decoding,
            false);
    // solve uv = v
    prpack_result* ret_v = solve_via_schur_gs(
            alpha,
            tol,
            num_vs,
            num_no_in_vs,
            num_no_out_vs,
            num_es,
            heads,
            tails,
            vals,
            ii,
            d,
            num_outlinks,
            v,
            encoding,
            decoding,
            false);
    // combine the u and v cases
    return combine_uv(num_vs, d, num_outlinks, encoding, alpha, ret_u, ret_v);
}

/** Gauss-Seidel using strongly connected components.
 * Notes:
 *   If not weighted, then we store x[i] = "x[i]/outdegree" to 
 *   avoid additional arithmetic.  We don't do this for the weighted
 *   case because the adjustment may not be constant.
 */
prpack_result* prpack_solver::solve_via_scc_gs(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es_inside,
        const int* heads_inside,
        const int* tails_inside,
        const double* vals_inside,
        const int num_es_outside,
        const int* heads_outside,
        const int* tails_outside,
        const double* vals_outside,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* uv,
        const int num_comps,
        const int* divisions,
        const int* encoding,
        const int* decoding,
        const bool should_normalize) {
    prpack_result* ret = new prpack_result();
    const bool weighted = vals_inside != NULL;
    // initialize uv values
    const double uv_const = 1.0/num_vs;
    const int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? prpack_utils::permute(num_vs, uv, encoding) : &uv_const;
    // CHECK initialize the solution with one iteration of GS from x=0.
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        x[i] = uv[uv_exists*i]/(1 - alpha*ii[i])/((weighted) ? 1 : num_outlinks[i]);
    // create x_outside
    double* x_outside = new double[num_vs];
    // run Gauss-Seidel for (I - alpha*P)*x = uv
    ret->num_es_touched = 0;
    for (int comp_i = 0; comp_i < num_comps; ++comp_i) {
        const int start_comp = divisions[comp_i];
        const int end_comp = (comp_i + 1 != num_comps) ? divisions[comp_i + 1] : num_vs;
        const bool parallelize = end_comp - start_comp > 512;
        // initialize relevant x_outside values
        for (int i = start_comp; i < end_comp; ++i) {
            x_outside[i] = 0;
            const int start_j = tails_outside[i];
            const int end_j = (i + 1 != num_vs) ? tails_outside[i + 1] : num_es_outside;
            for (int j = start_j; j < end_j; ++j)
                x_outside[i] += x[heads_outside[j]]*((weighted) ? vals_outside[j] : 1.);
            ret->num_es_touched += end_j - start_j;
        }
        double err, c;
        do {
            int num_es_touched = 0;
            err = c = 0;
            if (parallelize) {
                // iterate through vertices
                #pragma omp parallel for firstprivate(c) reduction(+:err, num_es_touched) schedule(dynamic, 64)
                for (int i = start_comp; i < end_comp; ++i) {
                    double new_val = x_outside[i];
                    const int start_j = tails_inside[i];
                    const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
                    if (weighted) {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]]*vals_inside[j];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    } else {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]*num_outlinks[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i])/num_outlinks[i];
                    }
                    num_es_touched += end_j - start_j;
                }
            } else {
                for (int i = start_comp; i < end_comp; ++i) {
                    double new_val = x_outside[i];
                    const int start_j = tails_inside[i];
                    const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
                    if (weighted) {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]]*vals_inside[j];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    } else {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]*num_outlinks[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i])/num_outlinks[i];
                    }
                    num_es_touched += end_j - start_j;
                }
            }
            // update iteration index
            ret->num_es_touched += num_es_touched;
        } while (err/(1 - alpha) >= tol*(end_comp - start_comp)/num_vs);
    }
    // undo num_outlinks transformation
    if (!weighted)
        for (int i = 0; i < num_vs; ++i)
            x[i] *= num_outlinks[i];
    // normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
    if (should_normalize)
        normalize(num_vs, x);
    // return results
    ret->x = prpack_utils::permute(num_vs, x, decoding);
    delete[] x;
    delete[] x_outside;
    if (uv_exists)
        delete[] uv;
    return ret;
}

prpack_result* prpack_solver::solve_via_scc_gs_uv(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es_inside,
        const int* heads_inside,
        const int* tails_inside,
        const double* vals_inside,
        const int num_es_outside,
        const int* heads_outside,
        const int* tails_outside,
        const double* vals_outside,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* u,
        const double* v,
        const int num_comps,
        const int* divisions,
        const int* encoding,
        const int* decoding) {
    // solve uv = u
    prpack_result* ret_u = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            vals_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            vals_outside,
            ii,
            d,
            num_outlinks,
            u,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // solve uv = v
    prpack_result* ret_v = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            vals_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            vals_outside,
            ii,
            d,
            num_outlinks,
            v,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // combine u and v
    return combine_uv(num_vs, d, num_outlinks, encoding, alpha, ret_u, ret_v);
}

// VARIOUS HELPER METHODS /////////////////////////////////////////////////////////////////////////

// Run Gaussian-Elimination (note: this changes A and returns the solution in b)
void prpack_solver::ge(const int sz, double* A, double* b) {
    // put into triangular form
    for (int i = 0, isz = 0; i < sz; ++i, isz += sz)
        for (int k = 0, ksz = 0; k < i; ++k, ksz += sz)
            if (A[isz + k] != 0) {
                const double coeff = A[isz + k]/A[ksz + k];
                A[isz + k] = 0;
                for (int j = k + 1; j < sz; ++j)
                    A[isz + j] -= coeff*A[ksz + j];
                b[i] -= coeff*b[k];
            }
    // backwards substitution
    for (int i = sz - 1, isz = (sz - 1)*sz; i >= 0; --i, isz -= sz) {
        for (int j = i + 1; j < sz; ++j)
            b[i] -= A[isz + j]*b[j];
        b[i] /= A[isz + i];
    }
}

// Normalize a vector to sum to 1.
void prpack_solver::normalize(const int length, double* x) {
    double norm = 0, c = 0;
    for (int i = 0; i < length; ++i) {
        COMPENSATED_SUM(norm, x[i], c);
    }
    norm = 1/norm;
    for (int i = 0; i < length; ++i)
        x[i] *= norm;
}

// Combine u and v results.
prpack_result* prpack_solver::combine_uv(
        const int num_vs,
        const double* d,
        const double* num_outlinks,
        const int* encoding,
        const double alpha,
        const prpack_result* ret_u,
        const prpack_result* ret_v) {
    prpack_result* ret = new prpack_result();
    const bool weighted = d != NULL;
    double delta_u = 0;
    double delta_v = 0;
    for (int i = 0; i < num_vs; ++i) {
        if ((weighted) ? (d[encoding[i]] == 1) : (num_outlinks[encoding[i]] < 0)) {
            delta_u += ret_u->x[i];
            delta_v += ret_v->x[i];
        }
    }
    const double s = ((1 - alpha)*alpha*delta_v)/(1 - alpha*delta_u);
    const double t = 1 - alpha;
    ret->x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        ret->x[i] = s*ret_u->x[i] + t*ret_v->x[i];
    ret->num_es_touched = ret_u->num_es_touched + ret_v->num_es_touched;
    // clean up and return
    delete ret_u;
    delete ret_v;
    return ret;
}

