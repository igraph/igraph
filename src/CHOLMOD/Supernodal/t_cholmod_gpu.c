/* ========================================================================== */
/* === Supernodal/t_cholmod_gpu ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Supernodal Module.  Copyright (C) 2005-2012, Timothy A. Davis
 * The CHOLMOD/Supernodal Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* GPU BLAS template routine for cholmod_super_numeric. */

/* ========================================================================== */
/* === include files and definitions ======================================== */
/* ========================================================================== */

#include "cholmod_template.h"

#undef L_ENTRY
#ifdef REAL
#define L_ENTRY 1
#else
#define L_ENTRY 2
#endif

/*
#define GPU_Printf  printf
*/
#define GPU_Printf

#define PAGE_SIZE (4*1024)
#define OK(cuda_operation) ((cuda_operation) == cudaSuccess)

/* ========================================================================== */
/* === gpu_init ============================================================= */
/* ========================================================================== */

void TEMPLATE (CHOLMOD (gpu_init))
(
    void *Cwork,
    Int maxSize,
    cholmod_common *Common
)
{
    Int i ;
    cublasStatus_t cublasError ;
    cudaError_t cudaErr ;
    size_t maxBytesSize, HostPinnedSize ;

    Common->GemmUsed = 0 ;

    GPU_Printf ("gpu_init : %p\n", (void *) ((size_t) Cwork & ~(PAGE_SIZE-1))) ;

    if (!(Common->cublasHandle))
    {

        /* ------------------------------------------------------------------ */
        /* create the CUDA BLAS handle */
        /* ------------------------------------------------------------------ */

        cublasError = cublasCreate (&(Common->cublasHandle)) ;
        if (cublasError != CUBLAS_STATUS_SUCCESS)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUBLAS initialization") ;
            return ;
        }

        /* ------------------------------------------------------------------ */
        /* create each CUDA stream */
        /* ------------------------------------------------------------------ */

        cudaErr = cudaStreamCreate (&(Common->cudaStreamSyrk)) ;
        if (cudaErr != cudaSuccess)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA stream initialization") ;
            return ;
        }

        cudaErr = cudaStreamCreate (&(Common->cudaStreamGemm)) ;
        if (cudaErr != cudaSuccess)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA stream initialization") ;
            return ;
        }

        cudaErr = cudaStreamCreate (&(Common->cudaStreamTrsm)) ;
        if (cudaErr != cudaSuccess)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA stream initialization") ;
            return ;
        }

        for (i = 0 ; i < 3 ; i++)
        {
            cudaErr = cudaStreamCreate (&(Common->cudaStreamPotrf [i])) ;
            if (cudaErr != cudaSuccess)
            {
                ERROR (CHOLMOD_GPU_PROBLEM, "CUDA stream initialization") ;
                return ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* create each CUDA event */
        /* ------------------------------------------------------------------ */

        for (i = 0 ; i < 2 ; i++)
        {
            cudaErr = cudaEventCreateWithFlags
                (&(Common->cublasEventPotrf [i]), cudaEventDisableTiming) ;
            if (cudaErr != cudaSuccess)
            {
                ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event") ;
                return ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* pin the Host memory */
    /* ---------------------------------------------------------------------- */

    Common->HostPinnedMemory = (void *) ((size_t) Cwork & ~(PAGE_SIZE-1)) ;
    maxBytesSize = sizeof (double)*L_ENTRY*maxSize ;

    /* Align on a 4K page boundary (it is no more necessary in 4.1 */
    HostPinnedSize  = 
        (((size_t) Cwork + maxBytesSize + PAGE_SIZE-1) & ~(PAGE_SIZE-1))
        - (size_t) (Common->HostPinnedMemory) ;

    GPU_Printf ("gpu HostPinnedSize: %g %p\n", (double) HostPinnedSize,
        Common->HostPinnedMemory) ;
    cudaErr = cudaHostRegister (Common->HostPinnedMemory,
        HostPinnedSize, 0) ;

    if (cudaErr != cudaSuccess)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUDA Pinning Memory") ;
        Common->HostPinnedMemory = NULL ;
    }
}


/* ========================================================================== */
/* === gpu_end ============================================================== */
/* ========================================================================== */

void TEMPLATE (CHOLMOD (gpu_end))
(
    cholmod_common *Common
)
{
    int i;
    /* unpin the Host memory */
    GPU_Printf ("gpu_end %p\n", Common->HostPinnedMemory) ;
    cudaError_t cudaErr = cudaHostUnregister (Common->HostPinnedMemory) ;
    if (cudaErr != cudaSuccess)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUDA Unpinning Memory") ;
        Common->HostPinnedMemory = NULL ;
    }
    /* ------------------------------------------------------------------ */
    /* destroy Cublas Handle */
    /* ------------------------------------------------------------------ */
    
    if (Common->cublasHandle) {
        cublasDestroy(Common->cublasHandle);
        Common->cublasHandle = NULL ;
    }
    /* ------------------------------------------------------------------ */
    /* destroy each CUDA stream */
    /* ------------------------------------------------------------------ */
    if (Common->cudaStreamSyrk) 
    {
        cudaStreamDestroy (Common->cudaStreamSyrk) ;
        Common->cudaStreamSyrk = NULL ;
    }
    if (Common->cudaStreamGemm) 
    {
        cudaStreamDestroy (Common->cudaStreamGemm) ;
    }
    if (Common->cudaStreamTrsm) 
    {
        cudaStreamDestroy (Common->cudaStreamTrsm) ;
        Common->cudaStreamTrsm = NULL ;
    }        

    for (i = 0 ; i < 3 ; i++)
    {
        if (Common->cudaStreamPotrf [i]) 
        {
            cudaStreamDestroy(Common->cudaStreamPotrf [i]) ;
            Common->cudaStreamPotrf [i] = NULL ;
        }
    }

    /* ------------------------------------------------------------------ */
    /* destroy each CUDA event */
    /* ------------------------------------------------------------------ */

    for (i = 0 ; i < 2 ; i++)
    {
        if (Common->cublasEventPotrf [i]) 
        {
            cudaEventDestroy( Common->cublasEventPotrf [i] ) ;
            Common->cublasEventPotrf [i] = NULL ;
        }
    }    
}


/* ========================================================================== */
/* === gpu_updateC ========================================================== */
/* ========================================================================== */

/* C = L (k1:n-1, kd1:kd2-1) * L (k1:k2-1, kd1:kd2-1)', except that k1:n-1
 * refers to all of the rows in L, but many of the rows are all zero.
 * Supernode d holds columns kd1 to kd2-1 of L.  Nonzero rows in the range
 * k1:k2-1 are in the list Ls [pdi1 ... pdi2-1], of size ndrow1.  Nonzero rows
 * in the range k2:n-1 are in the list Ls [pdi2 ... pdend], of size ndrow2.
 * Let L1 = L (Ls [pdi1 ... pdi2-1], kd1:kd2-1), and let L2 = L (Ls [pdi2 ...
 * pdend],  kd1:kd2-1).  C is ndrow2-by-ndrow1.  Let C1 be the first ndrow1
 * rows of C and let C2 be the last ndrow2-ndrow1 rows of C.  Only the lower
 * triangular part of C1 needs to be computed since C1 is symmetric.
 */

int TEMPLATE (CHOLMOD (gpu_updateC))
(
    Int ndrow1,         /* C is ndrow2-by-ndrow2 */   
    Int ndrow2,
    Int ndrow,          /* leading dimension of Lx */
    Int ndcol,          /* L1 is ndrow1-by-ndcol */
    Int pdx1,           /* L1 starts at Lx + L_ENTRY*pdx1 */
                        /* L2 starts at Lx + L_ENTRY*(pdx1 + ndrow1) */
    double *Lx,
    double *C,
    cholmod_common *Common
)
{
    double *devPtrLx, *devPtrC ;
    double alpha, beta ;
    cublasStatus_t cublasStatus ;
    cudaError_t cudaStat [2] ;
    Int ndrow3 ;

    Common->SyrkUsed = 0 ;
    Common->GemmUsed = 0 ;

    if ((ndrow2 < 512) || (ndcol <  128))
    {
        /* too small for the CUDA BLAS; use the CPU instead */
        return (0) ;
    }

    ndrow3 = ndrow2 - ndrow1 ;

#ifndef NTIMER
    Common->syrkStart = SuiteSparse_time ( ) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* allocate workspace on the GPU */
    /* ---------------------------------------------------------------------- */

    cudaStat [0] = cudaMalloc ((void **) &devPtrLx,
        ndrow2 * ndcol  * L_ENTRY * sizeof (devPtrLx [0])) ;
    cudaStat [1] = cudaMalloc ((void **) &devPtrC,
        ndrow2 * ndrow1 * L_ENTRY * sizeof (devPtrC [0])) ;
    Common->devSyrkGemmPtrLx = devPtrLx ;
    Common->devSyrkGemmPtrC  = devPtrC ;
    
    if (cudaStat [0] || cudaStat [1])
    {
        /* one or both cudaMalloc's failed */
        if (devPtrLx) cudaFree (devPtrLx) ;
        if (devPtrC)  cudaFree (devPtrC) ;
        GPU_Printf ("gpu malloc failed =%d,%d ndrow1=%d ndrow2=%d ndcol=%d\n",
            cudaStat [0], cudaStat [1], (int) ndrow1,
            (int) ndrow2, (int) ndcol) ;
        /* cudaMalloc failure is not an error, just bypass the GPU */
        return (0) ;
    }
    Common->SyrkUsed = 1 ;
#ifndef NTIMER
    Common->CHOLMOD_GPU_SYRK_CALLS++ ;
#endif

    /* ---------------------------------------------------------------------- */
    /* copy Lx to the GPU */
    /* ---------------------------------------------------------------------- */

    /* copy Lx in two steps on different streams.
     * (ldLx is shortened from ndrow to ndrow2) */
    cudaStat [0] = cudaMemcpy2DAsync (devPtrLx,
        ndrow2 * L_ENTRY * sizeof (devPtrLx [0]),
        Lx + L_ENTRY * pdx1, ndrow * L_ENTRY * sizeof (Lx [0]),
        ndrow1 * L_ENTRY * sizeof (devPtrLx [0]),
        ndcol, cudaMemcpyHostToDevice, Common->cudaStreamSyrk) ;
    if (cudaStat [0])
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
    }

    if (ndrow3 > 0)
    {
        Common->GemmUsed = 1 ;
        cudaStat [1] = cudaMemcpy2DAsync (devPtrLx + L_ENTRY*ndrow1,
            ndrow2 * L_ENTRY * sizeof (devPtrLx [0]),
            Lx + L_ENTRY * (pdx1 + ndrow1), ndrow * L_ENTRY * sizeof (Lx [0]),
            ndrow3 * L_ENTRY * sizeof (devPtrLx [0]),
            ndcol, cudaMemcpyHostToDevice, Common->cudaStreamGemm) ;
        if (cudaStat [1])
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* do the CUDA SYRK */
    /* ---------------------------------------------------------------------- */

    cublasStatus = cublasSetStream (Common->cublasHandle,
        Common->cudaStreamSyrk) ;
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
    }

    alpha  = 1.0 ;
    beta   = 0.0 ;
#ifdef REAL
    cublasStatus = cublasDsyrk (Common->cublasHandle,
        CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
        (int) ndrow1, (int) ndcol,          /* N, K: L1 is ndrow1-by-ndcol */
        &alpha,                             /* ALPHA:  1 */
        devPtrLx, ndrow2,                   /* A, LDA: L1, ndrow2 */
        &beta,                              /* BETA:   0 */
        devPtrC, ndrow2) ;                  /* C, LDC: C1 */
#else
    cublasStatus = cublasZherk (Common->cublasHandle,
        CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
        (int) ndrow1, (int) ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
        &alpha,                                 /* ALPHA:  1 */
        (const cuDoubleComplex *) devPtrLx, ndrow2, /* A, LDA: L1, ndrow2 */
        &beta,                                  /* BETA:   0 */
        (cuDoubleComplex *) devPtrC, ndrow2) ;  /* C, LDC: C1 */
#endif

    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
    }

    /* ---------------------------------------------------------------------- */
    /* partial copy of C to the GPU */
    /* ---------------------------------------------------------------------- */

    cudaStat [0] = cudaMemcpy2DAsync (C, ndrow2 * L_ENTRY * sizeof (C [0]),
        devPtrC, ndrow2 * L_ENTRY * sizeof (devPtrC [0]),
        ndrow1 * L_ENTRY * sizeof (devPtrC [0]),
        ndrow1, cudaMemcpyDeviceToHost, Common->cudaStreamSyrk) ;
    if (cudaStat [0])
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute remaining (ndrow2-ndrow1)-by-ndrow1 block of C, C2 = L2*L1' */
    /* ---------------------------------------------------------------------- */

    if (ndrow3 > 0)
    {
#ifndef REAL
        cuDoubleComplex calpha  = {1.0,0.0} ;
        cuDoubleComplex cbeta   = {0.0,0.0} ;
#endif

#ifndef NTIMER
        Common->CHOLMOD_GPU_GEMM_CALLS++ ;
#endif
        cublasStatus = cublasSetStream (Common->cublasHandle,
            Common->cudaStreamGemm) ;
        if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
        }

        /* ------------------------------------------------------------------ */
        /* do the CUDA BLAS dgemm */
        /* ------------------------------------------------------------------ */

#ifdef REAL
        alpha  = 1.0 ;
        beta   = 0.0 ;
        cublasStatus = cublasDgemm (Common->cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_T,
            ndrow3, ndrow1, ndcol,                  /* M, N, K */
            &alpha,                                 /* ALPHA:  1 */
            devPtrLx + L_ENTRY*(ndrow1),            /* A, LDA: L2, ndrow */
            ndrow2,
            devPtrLx,                               /* B, LDB: L1, ndrow */
            ndrow2,
            &beta,                                  /* BETA:   0 */
            devPtrC + L_ENTRY*ndrow1,               /* C, LDC: C2 */
            ndrow2) ;
#else
        cublasStatus = cublasZgemm (Common->cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_C,
            ndrow3, ndrow1, ndcol,                  /* M, N, K */
            &calpha,                                /* ALPHA:  1 */
            (const cuDoubleComplex *) devPtrLx + ndrow1, /* A, LDA: L2, ndrow */
            ndrow2,
            (const cuDoubleComplex *) devPtrLx,     /* B, LDB: L1, ndrow */
            ndrow2,
            &cbeta,                                 /* BETA:   0 */
            (cuDoubleComplex *)devPtrC + ndrow1,    /* C, LDC: C2 */
            ndrow2) ;
#endif

        if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
        }

        /* ------------------------------------------------------------------ */
        /* finish copy of C */
        /* ------------------------------------------------------------------ */

        cudaStat [0] = cudaMemcpy2DAsync (C + L_ENTRY*ndrow1,
            ndrow2 * L_ENTRY * sizeof (C [0]),
            devPtrC+ L_ENTRY*ndrow1, ndrow2 * L_ENTRY * sizeof (devPtrC [0]),
            ndrow3 * L_ENTRY * sizeof (devPtrC [0]),
            ndrow1, cudaMemcpyDeviceToHost, Common->cudaStreamGemm) ;
        if (cudaStat [0])
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
        }
    }

    return (1) ;
}


/* ========================================================================== */
/* === gpu_syncSyrk ========================================================= */
/* ========================================================================== */

/* synchronize with the CUDA BLAS dsyrk stream */

void TEMPLATE (CHOLMOD (gpu_syncSyrk))
(
    cholmod_common *Common
)
{
    if (Common->SyrkUsed)
    {
        cudaStreamSynchronize (Common->cudaStreamSyrk) ;
        if (!Common->GemmUsed)
        {
            cudaFree (Common->devSyrkGemmPtrLx) ;
            cudaFree (Common->devSyrkGemmPtrC) ;
            Common->devSyrkGemmPtrLx = NULL ;
            Common->devSyrkGemmPtrC = NULL ;
#ifndef NTIMER
            /* this actually sums time spend on Syrk and Gemm */
            Common->CHOLMOD_GPU_SYRK_TIME +=
                SuiteSparse_time ( ) -  Common->syrkStart ;
#endif
        }
    }
}


/* ========================================================================== */
/* === gpu_syncGemm ========================================================= */
/* ========================================================================== */

/* synchronize with the CUDA BLAS dgemm stream */

void TEMPLATE (CHOLMOD (gpu_syncGemm))
(
    cholmod_common *Common
)
{
    if (Common->GemmUsed)
    {
        cudaStreamSynchronize (Common->cudaStreamGemm) ;
        cudaFree (Common->devSyrkGemmPtrLx) ;
        cudaFree (Common->devSyrkGemmPtrC) ;
        Common->devSyrkGemmPtrLx = NULL ;
        Common->devSyrkGemmPtrC = NULL ;
#ifndef NTIMER
        /* this actually sums time spend on Syrk and Gemm */
        Common->CHOLMOD_GPU_SYRK_TIME +=
            SuiteSparse_time ( ) - Common->syrkStart ;
#endif
    }
}


/* ========================================================================== */
/* === gpu_lower_potrf ====================================================== */
/* ========================================================================== */

/* Cholesky factorzation (dpotrf) of a matrix S, operating on the lower
 * triangular part only.   S is nscol2-by-nscol2 with leading dimension nsrow.
 *
 * S is the top part of the supernode (the lower triangular matrx).
 * This function also copies the bottom rectangular part of the supernode (B)
 * onto the GPU, in preparation for gpu_triangular_solve. 
 */

int TEMPLATE (CHOLMOD (gpu_lower_potrf))
(
    Int nscol2,     /* S is nscol2-by-nscol2 */
    Int nsrow,      /* leading dimension of S */
    Int psx,        /* S is located at Lx + L_Entry*psx */
    double *Lx,     /* contains S; overwritten with Cholesky factor */
    Int *info,      /* BLAS info return value */
    cholmod_common *Common
)
{
    double *devPtrA, *devPtrB, *A ;
    double alpha, beta ;
    cudaError_t cudaStat ;
    cublasStatus_t cublasStatus ;
    Int j, nsrow2, nb, n, gpu_lda, lda, gpu_ldb ;
    int ilda, ijb, iinfo ;
#ifndef NTIMER
    double tstart = SuiteSparse_time ( ) ;
#endif

    if (nscol2 < 256)
    {
        /* too small for the CUDA BLAS; use the CPU instead */
        return (0) ;
    }

    nsrow2 = nsrow - nscol2 ;

    /* ---------------------------------------------------------------------- */
    /* heuristic to get the block size depending of the problem size */
    /* ---------------------------------------------------------------------- */

    nb = 128 ;
    if (nscol2 > 4096) nb = 256 ;
    if (nscol2 > 8192) nb = 384 ;
    n  = nscol2 ;
    gpu_lda = ((nscol2+31)/32)*32 ;
    lda = nsrow ;
    A = Lx + L_ENTRY*psx ;

    /* ---------------------------------------------------------------------- */
    /* free the dpotrf workspace, if allocated */
    /* ---------------------------------------------------------------------- */

    if (Common->devPotrfWork)
    {
        cudaFree (Common->devPotrfWork) ;
        Common->devPotrfWork = NULL ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine the GPU leading dimension of B */
    /* ---------------------------------------------------------------------- */

    gpu_ldb = 0 ;
    if (nsrow2 > 0)
    {
        gpu_ldb = ((nsrow2+31)/32)*32 ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate device memory for the factorization and for potential solve */
    /* ---------------------------------------------------------------------- */

    cudaStat = cudaMalloc ((void **) &devPtrA,
        gpu_lda * (gpu_lda + gpu_ldb) * L_ENTRY * sizeof (devPtrA [0])) ;
    if (cudaStat)
    {
        GPU_Printf ("@@gpu_lower_potrf cudaMalloc failed =%d gpu_lda=%d\n",
            cudaStat, (int) (gpu_lda)) ;
        /* cudaMalloc failure not fatal, GPU bypassed */
        return (0) ;
    }
#ifndef NTIMER
    Common->CHOLMOD_GPU_POTRF_CALLS++ ;
#endif

    /* ---------------------------------------------------------------------- */
    /* remember where device memory is, to be used by triangular solve later */
    /* ---------------------------------------------------------------------- */

    Common->devPotrfWork = devPtrA ;
    devPtrB = devPtrA + gpu_lda * gpu_lda * L_ENTRY ;

    /* ---------------------------------------------------------------------- */
    /* copy B in advance, for gpu_triangular_solve */
    /* ---------------------------------------------------------------------- */

    if (nsrow2 > 0)
    {
        cudaStat = cudaMemcpy2DAsync (devPtrB,
            gpu_ldb * L_ENTRY * sizeof (devPtrB [0]),
            Lx + L_ENTRY * (psx + nscol2),
            nsrow * L_ENTRY * sizeof (Lx [0]),
            nsrow2 * L_ENTRY * sizeof (devPtrB [0]),
            nscol2, cudaMemcpyHostToDevice, Common->cudaStreamTrsm) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* block Cholesky factorization of S */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < n ; j += nb)
    {
        Int jb = nb < (n-j) ? nb : (n-j) ;

        /* ------------------------------------------------------------------ */
        /* copy jb columns starting at the diagonal to the GPU */
        /* ------------------------------------------------------------------ */

        cudaStat = cudaMemcpy2DAsync (devPtrA + (j + j*gpu_lda)*L_ENTRY,
            gpu_lda * L_ENTRY * sizeof (devPtrA [0]),
            A + L_ENTRY*(j + j*lda),
            lda * L_ENTRY * sizeof (A [0]),
            (n-j) * L_ENTRY * sizeof (devPtrA [0]),
            jb, cudaMemcpyHostToDevice, Common->cudaStreamPotrf [0]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }

        /* ------------------------------------------------------------------ */
        /* define the dpotrf stream */
        /* ------------------------------------------------------------------ */

        cublasStatus = cublasSetStream (Common->cublasHandle,
            Common->cudaStreamPotrf [0]) ;
        if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
        }

        /* ------------------------------------------------------------------ */
        /* record the end of the copy of block L22 | L32 */
        /* ------------------------------------------------------------------ */

        cudaStat = cudaEventRecord (Common->cublasEventPotrf [0],
            Common->cudaStreamPotrf [0]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
        }

        /* ------------------------------------------------------------------ */
        /* do the CUDA BLAS dsyrk */
        /* ------------------------------------------------------------------ */

        alpha = -1.0 ;
        beta  = 1.0 ;
#ifdef REAL
        cublasStatus = cublasDsyrk (Common->cublasHandle,
            CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, jb, j,
            &alpha, devPtrA + j, gpu_lda,
            &beta,  devPtrA + j + j*gpu_lda, gpu_lda) ;
#else
        cublasStatus = cublasZherk (Common->cublasHandle,
            CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, jb, j,
            &alpha, (cuDoubleComplex*)devPtrA + j, gpu_lda,
            &beta,  (cuDoubleComplex*)devPtrA + j + j*gpu_lda, gpu_lda) ;
#endif
        if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
        }

        /* ------------------------------------------------------------------ */

        cudaStat = cudaEventRecord (Common->cublasEventPotrf [1],
            Common->cudaStreamPotrf [0]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
        }

        cudaStat = cudaStreamWaitEvent (Common->cudaStreamPotrf [1],
            Common->cublasEventPotrf [1], 0) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
        }

        /* ------------------------------------------------------------------ */
        /* copy back the jb columns on two different streams */
        /* ------------------------------------------------------------------ */

        cudaStat = cudaMemcpy2DAsync (A + L_ENTRY*(j + j*lda),
            lda * L_ENTRY * sizeof (double),
            devPtrA + L_ENTRY*(j + j*gpu_lda),
            gpu_lda * L_ENTRY * sizeof (double),
            L_ENTRY * sizeof (double)*jb, jb,
            cudaMemcpyDeviceToHost, Common->cudaStreamPotrf [1]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
        }

        cudaStat = cudaMemcpy2DAsync (A + L_ENTRY*j,
            lda * L_ENTRY * sizeof (double),
            devPtrA + L_ENTRY*j,
            gpu_lda * L_ENTRY * sizeof (double),
            L_ENTRY * sizeof (double)*jb, j,
            cudaMemcpyDeviceToHost, Common->cudaStreamPotrf [0]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }

        /* ------------------------------------------------------------------ */
        /* do the CUDA BLAS dgemm */
        /* ------------------------------------------------------------------ */

        if ((j+jb) < n)
        {
#ifdef REAL
            alpha = -1.0 ;
            beta  = 1.0 ;
            cublasStatus = cublasDgemm (Common->cublasHandle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                (n-j-jb), jb, j,
                &alpha,
                devPtrA + (j+jb), gpu_lda,
                devPtrA + (j)  , gpu_lda,
                &beta,
                devPtrA + (j+jb + j*gpu_lda), gpu_lda) ;
#else
            cuDoubleComplex calpha = {-1.0,0.0} ;
            cuDoubleComplex cbeta  = { 1.0,0.0} ;
            cublasStatus = cublasZgemm (Common->cublasHandle,
                CUBLAS_OP_N, CUBLAS_OP_C,
                (n-j-jb), jb, j,
                &calpha,
                (cuDoubleComplex*)devPtrA + (j+jb), gpu_lda,
                (cuDoubleComplex*)devPtrA + (j)  , gpu_lda,
                &cbeta,
                (cuDoubleComplex*)devPtrA + (j+jb + j*gpu_lda), gpu_lda) ;
#endif
            if (cublasStatus != CUBLAS_STATUS_SUCCESS)
            {
                ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
            }
        }

        cudaStat = cudaStreamSynchronize (Common->cudaStreamPotrf [1]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }

        /* ------------------------------------------------------------------ */
        /* compute the Cholesky factorization of the jbxjb block on the CPU */
        /* ------------------------------------------------------------------ */

        ilda = (int) lda ;
        ijb  = jb ;
#ifdef REAL
        LAPACK_DPOTRF ("L", &ijb, A + L_ENTRY * (j + j*lda), &ilda, &iinfo) ;
#else
        LAPACK_ZPOTRF ("L", &ijb, A + L_ENTRY * (j + j*lda), &ilda, &iinfo) ;
#endif
        *info = iinfo ;

        if (*info != 0)
        {
            *info = *info + j ;
            break ;
        }

        /* ------------------------------------------------------------------ */
        /* copy the result back to the GPU */
        /* ------------------------------------------------------------------ */

        cudaStat = cudaMemcpy2DAsync (devPtrA + L_ENTRY*(j + j*gpu_lda),
            gpu_lda * L_ENTRY * sizeof (double),
            A + L_ENTRY * (j + j*lda),
            lda * L_ENTRY * sizeof (double),
            L_ENTRY * sizeof (double) * jb, jb,
            cudaMemcpyHostToDevice, Common->cudaStreamPotrf [0]) ;
        if (cudaStat)
        {
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
        }

        /* ------------------------------------------------------------------ */
        /* do the CUDA BLAS dtrsm */
        /* ------------------------------------------------------------------ */

        if ((j+jb) < n)
        {
#ifdef REAL
             alpha  = 1.0 ;
             cublasStatus = cublasDtrsm (Common->cublasHandle,
                 CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                 CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
                 (n-j-jb), jb,
                 &alpha,
                 devPtrA + (j + j*gpu_lda), gpu_lda,
                 devPtrA + (j+jb + j*gpu_lda), gpu_lda) ;
#else
             cuDoubleComplex calpha  = {1.0,0.0};
             cublasStatus = cublasZtrsm (Common->cublasHandle,
                 CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                 CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT,
                 (n-j-jb), jb,
                 &calpha,
                 (cuDoubleComplex *)devPtrA + (j + j*gpu_lda), gpu_lda,
                 (cuDoubleComplex *)devPtrA + (j+jb + j*gpu_lda), gpu_lda) ;
#endif
            if (cublasStatus != CUBLAS_STATUS_SUCCESS)
            {
                ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
            }
        }
    }

    if (nsrow2 <= 0)
    {
        /* No TRSM necessary */
        cudaFree (Common->devPotrfWork) ;
        Common->devPotrfWork = NULL ;
    }

#ifndef NTIMER
    Common->CHOLMOD_GPU_POTRF_TIME += SuiteSparse_time ( ) - tstart ;
#endif
    return (1) ;
}


/* ========================================================================== */
/* === gpu_triangular_solve ================================================= */
/* ========================================================================== */

/* The current supernode is columns k1 to k2-1 of L.  Let L1 be the diagonal
 * block (factorized by dpotrf/zpotrf above; rows/cols k1:k2-1), and L2 be rows
 * k2:n-1 and columns k1:k2-1 of L.  The triangular system to solve is L2*L1' =
 * S2, where S2 is overwritten with L2.  More precisely, L2 = S2 / L1' in
 * MATLAB notation.
 */

/* Version with pre-allocation in POTRF */

int TEMPLATE (CHOLMOD (gpu_triangular_solve))
(
    Int nsrow2,     /* L1 and S2 are nsrow2-by-nscol2 */
    Int nscol2,     /* L1 is nscol2-by-nscol2 */
    Int nsrow,      /* leading dimension of L1, L2, and S2 */
    Int psx,        /* L1 is at Lx+L_ENTRY*psx; L2 at Lx+L_ENTRY*(psx+nscol2)*/
    double *Lx,     /* holds L1, L2, and S2 */
    cholmod_common *Common
)
{
    double *devPtrA, *devPtrB ;
    cudaError_t cudaStat ;
    cublasStatus_t cublasStatus ;
    Int gpu_lda, gpu_ldb ;
#ifdef REAL
    double alpha  = 1.0 ;
#else
    cuDoubleComplex calpha  = {1.0,0.0} ;
#endif

    if (!Common->devPotrfWork)
    {
        /* no workspace for triangular solve */
        return (0) ;
    }

#ifndef NTIMER
    double tstart = SuiteSparse_time ( ) ;
    Common->CHOLMOD_GPU_TRSM_CALLS++ ;
#endif

    gpu_lda = ((nscol2+31)/32)*32 ;
    gpu_ldb = ((nsrow2+31)/32)*32 ;

    devPtrA = Common->devPotrfWork ;
    devPtrB = devPtrA + gpu_lda * gpu_lda * L_ENTRY ;

    /* ---------------------------------------------------------------------- */
    /* start the trsm stream */
    /* ---------------------------------------------------------------------- */

    cublasStatus = cublasSetStream (Common->cublasHandle,
        Common->cudaStreamTrsm) ;
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
    }

    /* ---------------------------------------------------------------------- */
    /* do the CUDA BLAS dtrsm */
    /* ---------------------------------------------------------------------- */

#ifdef REAL
    cublasStatus = cublasDtrsm (Common->cublasHandle,
        CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
        CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
        nsrow2, nscol2,                             /* M, N */
        &alpha,                                     /* ALPHA:  1 */
        devPtrA, gpu_lda,                           /* A, LDA */
        devPtrB, gpu_ldb) ;                         /* B, LDB */
#else
    cublasStatus = cublasZtrsm (Common->cublasHandle,
        CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
        CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT,
        nsrow2, nscol2,                             /* M, N */
        &calpha,                                    /* ALPHA:  1 */
        (const cuDoubleComplex *) devPtrA, gpu_lda, /* A, LDA */
        (cuDoubleComplex *) devPtrB, gpu_ldb) ;     /* B, LDB: nsrow2 */
#endif
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
    }

    /* ---------------------------------------------------------------------- */
    /* copy result back to the CPU */
    /* ---------------------------------------------------------------------- */

    cudaStat = cudaMemcpy2DAsync (Lx + L_ENTRY*(psx + nscol2),
        nsrow * L_ENTRY * sizeof (Lx [0]),
        devPtrB, gpu_ldb * L_ENTRY * sizeof (devPtrB [0]),
        nsrow2 * L_ENTRY * sizeof (devPtrB [0]),
        nscol2, cudaMemcpyDeviceToHost, Common->cudaStreamTrsm) ;
    if (cudaStat)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
    }

    /* ---------------------------------------------------------------------- */
    /* synchronize with the GPU */
    /* ---------------------------------------------------------------------- */

    cudaStat = cudaThreadSynchronize ( ) ;
    if (cudaStat)
    {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU synchronization failure") ;
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace and return */
    /* ---------------------------------------------------------------------- */
    
    cudaFree (Common->devPotrfWork) ;
    Common->devPotrfWork = NULL ;
#ifndef NTIMER
    Common->CHOLMOD_GPU_TRSM_TIME += SuiteSparse_time ( ) - tstart ;
#endif
    return (1) ;
}

#undef REAL
#undef COMPLEX
#undef ZOMPLEX
