#ifndef PRPACK_SOLVER
#define PRPACK_SOLVER
#include "prpack_base_graph.h"
#include "prpack_csc.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include "prpack_preprocessed_ge_graph.h"
#include "prpack_preprocessed_gs_graph.h"
#include "prpack_preprocessed_scc_graph.h"
#include "prpack_preprocessed_schur_graph.h"
#include "prpack_result.h"

// TODO Make this a user configurable variable
#define PRPACK_SOLVER_MAX_ITERS 1000000

namespace prpack {

    // Solver class.
    class prpack_solver {
        private:
            // instance variables
            double read_time;
            prpack_base_graph* bg;
            prpack_preprocessed_ge_graph* geg;
            prpack_preprocessed_gs_graph* gsg;
            prpack_preprocessed_schur_graph* sg;
            prpack_preprocessed_scc_graph* sccg;
			bool owns_bg;
            // methods
            void initialize();
            static prpack_result* solve_via_ge(
                    const double alpha,
                    const double tol,
                    const int num_vs,
                    const double* matrix,
                    const double* uv);
            static prpack_result* solve_via_ge_uv(
                    const double alpha,
                    const double tol,
                    const int num_vs,
                    const double* matrix,
                    const double* d,
                    const double* u,
                    const double* v);
            static prpack_result* solve_via_gs(
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
                    const double* v);
            static prpack_result* solve_via_gs_err(
                    const double alpha,
                    const double tol,
                    const int num_vs,
                    const int num_es,
                    const int* heads,
                    const int* tails,
                    const double* ii,
                    const double* num_outlinks,
                    const double* u,
                    const double* v);
            static prpack_result* solve_via_schur_gs(
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
                    const bool should_normalize = true);
            static prpack_result* solve_via_schur_gs_uv(
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
                    const int* decoding);
            static prpack_result* solve_via_scc_gs(
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
                    const bool should_normalize = true);
            static prpack_result* solve_via_scc_gs_uv(
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
                    const int* decoding);
            static void ge(const int sz, double* A, double* b);
            static void normalize(const int length, double* x);
            static prpack_result* combine_uv(
                    const int num_vs,
                    const double* d,
                    const double* num_outlinks,
                    const int* encoding,
                    const double alpha,
                    const prpack_result* ret_u,
                    const prpack_result* ret_v);
        public:
            // constructors
            prpack_solver(const prpack_csc* g);
            prpack_solver(const prpack_int64_csc* g);
            prpack_solver(const prpack_csr* g);
            prpack_solver(const prpack_edge_list* g);
            prpack_solver(prpack_base_graph* g, bool owns_bg=true);
            prpack_solver(const char* filename, const char* format, const bool weighted);
            // destructor
            ~prpack_solver();
            // methods
            int get_num_vs();
            prpack_result* solve(const double alpha, const double tol, const char* method);
            prpack_result* solve(
                    const double alpha,
                    const double tol,
                    const double* u,
                    const double* v,
                    const char* method);
    };

}

#endif
