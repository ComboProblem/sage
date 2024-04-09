from libcpp cimport bool as cppbool
from libcpp.vector cimport vector as cppvector
from gmpy2 cimport mpz # so we hope this works.

# integers
cdef extern from "gmp.h":
    # gmp integer
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_limb_signed_t
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr
    void mpz_init(mpz_t)
    void mpz_clear(mpz_t)
    cdef mpz_t* address_of_mpz "&"(mpz_t x)

cdef extern from "gmpxx.h":
    # gmp integer
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)

# from the sage libary
ctypedef mp_limb_t ulong
ctypedef mp_limb_signed_t slong

cdef extern from "flint/fmpz.h":
    # flint to gmp
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]
    void fmpz_get_mpz(mpz_t x, const fmpz_t f) noexcept
    void fmpz_init(fmpz_t f)
    void fmpz_set_ui(fmpz_t f, ulong g)
    void fmpz_clear(fmpz_t f)
    int fmpz_print(const fmpz_t x)

cdef extern from "pplite/pplite.hh" namespace "pplite":
    pass

cdef extern from "pplite/globals.hh" namespace "pplite":
    ctypedef size_t dim_type

cdef extern from "pplite/Var.hh" namespace "pplite":
    cdef cppclass Var
    cdef cppclass Var:
        Var(dim_type i)
        dim_type id()
        dim_type space_dim()
    ctypedef struct Vars_Set

cdef extern from "pplite/FLINT_Integer.hh" namespace "pplite":
    cdef cppclass FLINT_Integer
    cdef cppclass FLINT_Integer:
        FLINT_Integer()
        FLINT_Integer(const fmpz_t z)
        void print()

    
cdef extern from "pplite/Linear_Expr.hh" namespace "pplite":
    ctypedef cppvector[FLINT_Integer] Impl
    cdef cppclass Linear_Expr
    cdef cppclass Linear_Expr:
        Linear_Expr()
        Linear_Expr(Var v)
        Linear_Expr(Linear_Expr &e)
        dim_type id()  #methods? #linear_expr.hh lines 39-112
        dim_type space_dim()
        void set_space_dim(dim_type dim)
        void swap_space_dims(dim_type i, dim_type j)
        void shift_space_dims(dim_type start, dim_type n)
        void shift_space_dims(Var start, dim_type n)
        Impl impl()


#        void permute_space_dims_cycle(const Dims& cycle, dim_type d)
        
#cdef extern from "poly.hh" namespace "pplite":
#    cdeftype struct Poly_Impl
#    cdef cppclass Poly

# cdef extern from "PolySet.hh" namespace "pplite":
#     cdef cppclass PloySet

# cdef extern from "Poly_Rel.hh" namespace "pplite":
#     cdef cppclass Poly_Con_Rel
#     cdef cppclass Poly_Gen_Rel

# cdef extern from "Gen.hh" namespace "pplite":
#     cdef cppclass Gen

# cdef extern from "Con.hh" namespace "pplite":
    # cdef cppclass Con

#do I need utils.hh
        

# cdef extern from  "Two_Poly.hh" namespace "pplite":
#     cdef cppclass Two_Poly


#sclar_prod.hh looks like it deals with scalar products, has helper functions for other classes and verifications if a point is inndeed in an NNC polyhedron

# rational_fwd.hh gmp vs flint rationals
# rational.hh usage of gmp vs flint
# poly_winden.hh helper funcitons for widing operator?
# poly_templ.hh this does something but what I'm not sure. i think it deals with lines and rays of a polyhedron. 
# poly_min.hh internal file for double description method defineds helper fucntions for DD method
# poly_stats.hh testing and methods to describes properies of NNC polyhedrons
# poly_rel.hh ?????
# polyset_templ.hh internal stuff as a templeate for polyset? Its polyssible 
# PolySet.hh seems to do stuff important, possibly not neeeded. it seems to be all wrapped int eh poly_impl struc for data

# cdef extern from "PolySet.hh" namespace "pplite":
#     cdef cppclass PolySet

# cdef extern from "Poly.hh" namespace "pplite":
#     cdef cppclass Poly
#     ctypedef struct Poly_Impl
    
# # low_level_stats.hh can ignore
# # local_stats.hh is for testing and timing purposus. 
# # bits.hh doesn't look like it has useful fuserfacing code
# # B_Poly.hh safe to ignore i think 



    
# cdef extern from "Con.hh" namespace "pplite":
#     cdef cppclass Con
    
# cdef extern from "Affine_Expr.hh" namespace "pplite":
#     cdef cppclass Affine_Expr 
    