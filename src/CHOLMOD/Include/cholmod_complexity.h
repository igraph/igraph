/* ========================================================================== */
/* === Include/cholmod_complexity.h ========================================= */
/* ========================================================================== */

/* Define operations on pattern, real, complex, and zomplex objects.
 *
 * The xtype of an object defines it numerical type.  A qttern object has no
 * numerical values (A->x and A->z are NULL).  A real object has no imaginary
 * qrt (A->x is used, A->z is NULL).  A complex object has an imaginary qrt
 * that is stored interleaved with its real qrt (A->x is of size 2*nz, A->z
 * is NULL).  A zomplex object has both real and imaginary qrts, which are
 * stored seqrately, as in MATLAB (A->x and A->z are both used).
 *
 * XTYPE is CHOLMOD_PATTERN, _REAL, _COMPLEX or _ZOMPLEX, and is the xtype of
 * the template routine under construction.  XTYPE2 is equal to XTYPE, except
 * if XTYPE is CHOLMOD_PATTERN, in which case XTYPE is CHOLMOD_REAL.
 * XTYPE and XTYPE2 are defined in cholmod_template.h.  
 */

/* -------------------------------------------------------------------------- */
/* pattern */
/* -------------------------------------------------------------------------- */

#define P_TEMPLATE(name)		p_ ## name
#define P_ASSIGN2(x,z,p,ax,az,q)	x [p] = 1
#define P_PRINT(k,x,z,p)		PRK(k, ("1"))

/* -------------------------------------------------------------------------- */
/* real */
/* -------------------------------------------------------------------------- */

#define R_TEMPLATE(name)			r_ ## name
#define R_ASSEMBLE(x,z,p,ax,az,q)		x [p] += ax [q]
#define R_ASSIGN(x,z,p,ax,az,q)			x [p]  = ax [q]
#define R_ASSIGN_CONJ(x,z,p,ax,az,q)		x [p]  = ax [q]
#define R_ASSIGN_REAL(x,p,ax,q)			x [p]  = ax [q]
#define R_XTYPE_OK(type)			((type) == CHOLMOD_REAL)
#define R_IS_NONZERO(ax,az,q)			IS_NONZERO (ax [q])
#define R_IS_ZERO(ax,az,q)			IS_ZERO (ax [q])
#define R_IS_ONE(ax,az,q)			(ax [q] == 1)
#define R_MULT(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] * bx [r]
#define R_MULTADD(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define R_MULTSUB(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define R_MULTADDCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define R_MULTSUBCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define R_ADD(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] + bx [r]
#define R_ADD_REAL(x,p, ax,q, bx,r)		x [p]  = ax [q] + bx [r]
#define R_CLEAR(x,z,p)				x [p]  = 0
#define R_CLEAR_IMAG(x,z,p)
#define R_DIV(x,z,p,ax,az,q)			x [p] /= ax [q]
#define R_LLDOT(x,p, ax,az,q)			x [p] -= ax [q] * ax [q]
#define R_PRINT(k,x,z,p)			PRK(k, ("%24.16e", x [p]))

#define R_DIV_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] / bx [r]
#define R_MULT_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] * bx [r]

#define R_LDLDOT(x,p, ax,az,q, bx,r)		x [p] -=(ax[q] * ax[q])/ bx[r]

/* -------------------------------------------------------------------------- */
/* complex */
/* -------------------------------------------------------------------------- */

#define C_TEMPLATE(name)		c_ ## name
#define CT_TEMPLATE(name)		ct_ ## name

#define C_ASSEMBLE(x,z,p,ax,az,q) \
    x [2*(p)  ] += ax [2*(q)  ] ; \
    x [2*(p)+1] += ax [2*(q)+1]

#define C_ASSIGN(x,z,p,ax,az,q) \
    x [2*(p)  ] = ax [2*(q)  ] ; \
    x [2*(p)+1] = ax [2*(q)+1]

#define C_ASSIGN_REAL(x,p,ax,q)			x [2*(p)]  = ax [2*(q)]

#define C_ASSIGN_CONJ(x,z,p,ax,az,q) \
    x [2*(p)  ] =  ax [2*(q)  ] ; \
    x [2*(p)+1] = -ax [2*(q)+1]

#define C_XTYPE_OK(type)		((type) == CHOLMOD_COMPLEX)

#define C_IS_NONZERO(ax,az,q) \
    (IS_NONZERO (ax [2*(q)]) || IS_NONZERO (ax [2*(q)+1]))

#define C_IS_ZERO(ax,az,q) \
    (IS_ZERO (ax [2*(q)]) && IS_ZERO (ax [2*(q)+1]))

#define C_IS_ONE(ax,az,q) \
    ((ax [2*(q)] == 1) && IS_ZERO (ax [2*(q)+1]))

#define C_IMAG_IS_NONZERO(ax,az,q)  (IS_NONZERO (ax [2*(q)+1]))

#define C_MULT(x,z,p, ax,az,q, bx,bz,r) \
x [2*(p)  ] = ax [2*(q)  ] * bx [2*(r)] - ax [2*(q)+1] * bx [2*(r)+1] ; \
x [2*(p)+1] = ax [2*(q)+1] * bx [2*(r)] + ax [2*(q)  ] * bx [2*(r)+1]

#define C_MULTADD(x,z,p, ax,az,q, bx,bz,r) \
x [2*(p)  ] += ax [2*(q)  ] * bx [2*(r)] - ax [2*(q)+1] * bx [2*(r)+1] ; \
x [2*(p)+1] += ax [2*(q)+1] * bx [2*(r)] + ax [2*(q)  ] * bx [2*(r)+1]

#define C_MULTSUB(x,z,p, ax,az,q, bx,bz,r) \
x [2*(p)  ] -= ax [2*(q)  ] * bx [2*(r)] - ax [2*(q)+1] * bx [2*(r)+1] ; \
x [2*(p)+1] -= ax [2*(q)+1] * bx [2*(r)] + ax [2*(q)  ] * bx [2*(r)+1]

/* s += conj(a)*b */
#define C_MULTADDCONJ(x,z,p, ax,az,q, bx,bz,r) \
x [2*(p)  ] +=   ax [2*(q)  ]  * bx [2*(r)] + ax [2*(q)+1] * bx [2*(r)+1] ; \
x [2*(p)+1] += (-ax [2*(q)+1]) * bx [2*(r)] + ax [2*(q)  ] * bx [2*(r)+1]

/* s -= conj(a)*b */
#define C_MULTSUBCONJ(x,z,p, ax,az,q, bx,bz,r) \
x [2*(p)  ] -=   ax [2*(q)  ]  * bx [2*(r)] + ax [2*(q)+1] * bx [2*(r)+1] ; \
x [2*(p)+1] -= (-ax [2*(q)+1]) * bx [2*(r)] + ax [2*(q)  ] * bx [2*(r)+1]

#define C_ADD(x,z,p, ax,az,q, bx,bz,r) \
    x [2*(p)  ] = ax [2*(q)  ] + bx [2*(r)  ] ; \
    x [2*(p)+1] = ax [2*(q)+1] + bx [2*(r)+1]

#define C_ADD_REAL(x,p, ax,q, bx,r) \
    x [2*(p)] = ax [2*(q)] + bx [2*(r)]

#define C_CLEAR(x,z,p) \
    x [2*(p)  ] = 0 ; \
    x [2*(p)+1] = 0

#define C_CLEAR_IMAG(x,z,p) \
    x [2*(p)+1] = 0

/* s = s / a */
#define C_DIV(x,z,p,ax,az,q) \
    Common->complex_divide ( \
	      x [2*(p)],  x [2*(p)+1], \
	     ax [2*(q)], ax [2*(q)+1], \
	     &x [2*(p)], &x [2*(p)+1])

/* s -= conj(a)*a ; note that the result of conj(a)*a is real */
#define C_LLDOT(x,p, ax,az,q) \
    x [2*(p)] -= ax [2*(q)] * ax [2*(q)] + ax [2*(q)+1] * ax [2*(q)+1]

#define C_PRINT(k,x,z,p) PRK(k, ("(%24.16e,%24.16e)", x [2*(p)], x [2*(p)+1]))

#define C_DIV_REAL(x,z,p, ax,az,q, bx,r) \
    x [2*(p)  ] = ax [2*(q)  ] / bx [2*(r)] ; \
    x [2*(p)+1] = ax [2*(q)+1] / bx [2*(r)]

#define C_MULT_REAL(x,z,p, ax,az,q, bx,r) \
    x [2*(p)  ] = ax [2*(q)  ] * bx [2*(r)] ; \
    x [2*(p)+1] = ax [2*(q)+1] * bx [2*(r)]

/* s -= conj(a)*a/t */
#define C_LDLDOT(x,p, ax,az,q, bx,r) \
    x [2*(p)] -= (ax [2*(q)] * ax [2*(q)] + ax [2*(q)+1] * ax [2*(q)+1]) / bx[r]

/* -------------------------------------------------------------------------- */
/* zomplex */
/* -------------------------------------------------------------------------- */

#define Z_TEMPLATE(name)		z_ ## name
#define ZT_TEMPLATE(name)		zt_ ## name

#define Z_ASSEMBLE(x,z,p,ax,az,q) \
    x [p] += ax [q] ; \
    z [p] += az [q]

#define Z_ASSIGN(x,z,p,ax,az,q) \
    x [p] = ax [q] ; \
    z [p] = az [q]

#define Z_ASSIGN_REAL(x,p,ax,q)			x [p]  = ax [q]

#define Z_ASSIGN_CONJ(x,z,p,ax,az,q) \
    x [p] =  ax [q] ; \
    z [p] = -az [q]

#define Z_XTYPE_OK(type)		((type) == CHOLMOD_ZOMPLEX)

#define Z_IS_NONZERO(ax,az,q) \
    (IS_NONZERO (ax [q]) || IS_NONZERO (az [q]))

#define Z_IS_ZERO(ax,az,q) \
    (IS_ZERO (ax [q]) && IS_ZERO (az [q]))

#define Z_IS_ONE(ax,az,q) \
    ((ax [q] == 1) && IS_ZERO (az [q]))

#define Z_IMAG_IS_NONZERO(ax,az,q)  (IS_NONZERO (az [q]))

#define Z_MULT(x,z,p, ax,az,q, bx,bz,r) \
    x [p] = ax [q] * bx [r] - az [q] * bz [r] ; \
    z [p] = az [q] * bx [r] + ax [q] * bz [r]

#define Z_MULTADD(x,z,p, ax,az,q, bx,bz,r) \
    x [p] += ax [q] * bx [r] - az [q] * bz [r] ; \
    z [p] += az [q] * bx [r] + ax [q] * bz [r]

#define Z_MULTSUB(x,z,p, ax,az,q, bx,bz,r) \
    x [p] -= ax [q] * bx [r] - az [q] * bz [r] ; \
    z [p] -= az [q] * bx [r] + ax [q] * bz [r]

#define Z_MULTADDCONJ(x,z,p, ax,az,q, bx,bz,r) \
    x [p] +=   ax [q]  * bx [r] + az [q] * bz [r] ; \
    z [p] += (-az [q]) * bx [r] + ax [q] * bz [r]

#define Z_MULTSUBCONJ(x,z,p, ax,az,q, bx,bz,r) \
    x [p] -=   ax [q]  * bx [r] + az [q] * bz [r] ; \
    z [p] -= (-az [q]) * bx [r] + ax [q] * bz [r]

#define Z_ADD(x,z,p, ax,az,q, bx,bz,r) \
	x [p] = ax [q] + bx [r] ; \
	z [p] = az [q] + bz [r]

#define Z_ADD_REAL(x,p, ax,q, bx,r) \
	x [p] = ax [q] + bx [r]

#define Z_CLEAR(x,z,p) \
    x [p] = 0 ; \
    z [p] = 0

#define Z_CLEAR_IMAG(x,z,p) \
    z [p] = 0

/* s = s/a */
#define Z_DIV(x,z,p,ax,az,q) \
    Common->complex_divide (x [p], z [p], ax [q], az [q], &x [p], &z [p])

/* s -= conj(a)*a ; note that the result of conj(a)*a is real */
#define Z_LLDOT(x,p, ax,az,q) \
    x [p] -= ax [q] * ax [q] + az [q] * az [q]

#define Z_PRINT(k,x,z,p)	PRK(k, ("(%24.16e,%24.16e)", x [p], z [p]))

#define Z_DIV_REAL(x,z,p, ax,az,q, bx,r) \
    x [p] = ax [q] / bx [r] ; \
    z [p] = az [q] / bx [r]

#define Z_MULT_REAL(x,z,p, ax,az,q, bx,r) \
    x [p] = ax [q] * bx [r] ; \
    z [p] = az [q] * bx [r]

/* s -= conj(a)*a/t */
#define Z_LDLDOT(x,p, ax,az,q, bx,r) \
    x [p] -= (ax [q] * ax [q] + az [q] * az [q]) / bx[r]

/* -------------------------------------------------------------------------- */
/* all classes */
/* -------------------------------------------------------------------------- */

/* Check if A->xtype and the two arrays A->x and A->z are valid.  Set status to
 * invalid, unless status is already "out of memory".  A can be a sparse matrix,
 * dense matrix, factor, or triplet. */

#define RETURN_IF_XTYPE_INVALID(A,xtype1,xtype2,result) \
{ \
    if ((A)->xtype < (xtype1) || (A)->xtype > (xtype2) || \
        ((A)->xtype != CHOLMOD_PATTERN && ((A)->x) == NULL) || \
	((A)->xtype == CHOLMOD_ZOMPLEX && ((A)->z) == NULL)) \
    { \
	if (Common->status != CHOLMOD_OUT_OF_MEMORY) \
	{ \
	    ERROR (CHOLMOD_INVALID, "invalid xtype") ; \
	} \
	return (result) ; \
    } \
}
