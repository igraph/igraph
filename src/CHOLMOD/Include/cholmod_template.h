/* ========================================================================== */
/* === Include/cholmod_template.h =========================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* undefine current xtype macros, and then define macros for current type */
/* -------------------------------------------------------------------------- */

#undef TEMPLATE
#undef XTYPE
#undef XTYPE2
#undef XTYPE_OK
#undef ENTRY_IS_NONZERO
#undef ENTRY_IS_ZERO
#undef ENTRY_IS_ONE
#undef IMAG_IS_NONZERO

#undef ASSEMBLE
#undef ASSIGN
#undef ASSIGN_CONJ
#undef ASSIGN2
#undef ASSIGN2_CONJ
#undef ASSIGN_REAL
#undef MULT
#undef MULTADD
#undef ADD
#undef ADD_REAL
#undef MULTSUB
#undef MULTADDCONJ
#undef MULTSUBCONJ
#undef LLDOT
#undef CLEAR
#undef DIV
#undef DIV_REAL
#undef MULT_REAL
#undef CLEAR_IMAG
#undef LDLDOT
#undef PREFIX

#undef ENTRY_SIZE

#undef XPRINT0
#undef XPRINT1
#undef XPRINT2
#undef XPRINT3

/* -------------------------------------------------------------------------- */
/* pattern */
/* -------------------------------------------------------------------------- */


#ifdef PATTERN

#define PREFIX				    p_
#define TEMPLATE(name)			    P_TEMPLATE(name)
#define XTYPE				    CHOLMOD_PATTERN
#define XTYPE2				    CHOLMOD_REAL
#define XTYPE_OK(type)			    (TRUE)
#define ENTRY_IS_NONZERO(ax,az,q)	    (TRUE)
#define ENTRY_IS_ZERO(ax,az,q)		    (FALSE)
#define ENTRY_IS_ONE(ax,az,q)		    (TRUE)
#define IMAG_IS_NONZERO(ax,az,q)	    (FALSE)
#define ENTRY_SIZE			    0

#define ASSEMBLE(x,z,p,ax,az,q)
#define ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    P_ASSIGN2(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    P_ASSIGN2(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)
#define CLEAR_IMAG(x,z,p)
#define DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    P_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    P_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    P_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    P_PRINT(3,x,z,p)

/* -------------------------------------------------------------------------- */
/* real */
/* -------------------------------------------------------------------------- */

#elif defined (REAL)

#define PREFIX				    r_
#define TEMPLATE(name)			    R_TEMPLATE(name)
#define XTYPE				    CHOLMOD_REAL
#define XTYPE2				    CHOLMOD_REAL
#define XTYPE_OK(type)			    R_XTYPE_OK(type)
#define ENTRY_IS_NONZERO(ax,az,q)	    R_IS_NONZERO(ax,az,q)
#define ENTRY_IS_ZERO(ax,az,q)		    R_IS_ZERO(ax,az,q)
#define ENTRY_IS_ONE(ax,az,q)		    R_IS_ONE(ax,az,q)
#define IMAG_IS_NONZERO(ax,az,q)	    (FALSE)
#define ENTRY_SIZE			    1

#define ASSEMBLE(x,z,p,ax,az,q)		    R_ASSEMBLE(x,z,p,ax,az,q) 
#define ASSIGN(x,z,p,ax,az,q)		    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)	    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)		    R_ASSIGN_REAL(x,p,ax,q)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)	    R_MULT(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)     R_MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)	    R_ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)	    R_ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)     R_MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    R_MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    R_MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)		    R_LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)			    R_CLEAR(x,z,p) 
#define CLEAR_IMAG(x,z,p)		    R_CLEAR_IMAG(x,z,p) 
#define DIV(x,z,p,ax,az,q)		    R_DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)	    R_DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)	    R_MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)	    R_LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    R_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    R_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    R_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    R_PRINT(3,x,z,p)

/* -------------------------------------------------------------------------- */
/* complex */
/* -------------------------------------------------------------------------- */

#elif defined (COMPLEX)

#define PREFIX				    c_

#ifdef NCONJUGATE
#define TEMPLATE(name)			    CT_TEMPLATE(name)
#else
#define TEMPLATE(name)			    C_TEMPLATE(name)
#endif

#define ASSEMBLE(x,z,p,ax,az,q)		    C_ASSEMBLE(x,z,p,ax,az,q) 
#define ASSIGN(x,z,p,ax,az,q)		    C_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)	    C_ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    C_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    C_ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)		    C_ASSIGN_REAL(x,p,ax,q)
#define XTYPE				    CHOLMOD_COMPLEX
#define XTYPE2				    CHOLMOD_COMPLEX
#define XTYPE_OK(type)			    C_XTYPE_OK(type)
#define ENTRY_IS_NONZERO(ax,az,q)	    C_IS_NONZERO(ax,az,q)
#define ENTRY_IS_ZERO(ax,az,q)		    C_IS_ZERO(ax,az,q)
#define ENTRY_IS_ONE(ax,az,q)		    C_IS_ONE(ax,az,q)
#define IMAG_IS_NONZERO(ax,az,q)	    C_IMAG_IS_NONZERO(ax,az,q)
#define ENTRY_SIZE			    2

#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)     C_MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)	    C_MULT(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)	    C_ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)	    C_ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)     C_MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    C_MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    C_MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)		    C_LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)			    C_CLEAR(x,z,p) 
#define CLEAR_IMAG(x,z,p)		    C_CLEAR_IMAG(x,z,p) 
#define DIV(x,z,p,ax,az,q)		    C_DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)	    C_DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)	    C_MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)	    C_LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    C_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    C_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    C_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    C_PRINT(3,x,z,p)

/* -------------------------------------------------------------------------- */
/* zomplex */
/* -------------------------------------------------------------------------- */

#elif defined (ZOMPLEX)

#define PREFIX				    z_

#ifdef NCONJUGATE
#define TEMPLATE(name)			    ZT_TEMPLATE(name)
#else
#define TEMPLATE(name)			    Z_TEMPLATE(name)
#endif

#define ASSEMBLE(x,z,p,ax,az,q)		    Z_ASSEMBLE(x,z,p,ax,az,q) 
#define ASSIGN(x,z,p,ax,az,q)		    Z_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)	    Z_ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    Z_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    Z_ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)		    Z_ASSIGN_REAL(x,p,ax,q)
#define XTYPE				    CHOLMOD_ZOMPLEX
#define XTYPE2				    CHOLMOD_ZOMPLEX
#define XTYPE_OK(type)			    Z_XTYPE_OK(type)
#define ENTRY_IS_NONZERO(ax,az,q)	    Z_IS_NONZERO(ax,az,q)
#define ENTRY_IS_ZERO(ax,az,q)		    Z_IS_ZERO(ax,az,q)
#define ENTRY_IS_ONE(ax,az,q)		    Z_IS_ONE(ax,az,q)
#define IMAG_IS_NONZERO(ax,az,q)	    Z_IMAG_IS_NONZERO(ax,az,q)
#define ENTRY_SIZE			    1

#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)     Z_MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)	    Z_MULT(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)	    Z_ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)	    Z_ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)     Z_MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    Z_MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    Z_MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)		    Z_LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)			    Z_CLEAR(x,z,p) 
#define CLEAR_IMAG(x,z,p)		    Z_CLEAR_IMAG(x,z,p) 
#define DIV(x,z,p,ax,az,q)		    Z_DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)	    Z_DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)	    Z_MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)	    Z_LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    Z_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    Z_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    Z_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    Z_PRINT(3,x,z,p)

#endif
