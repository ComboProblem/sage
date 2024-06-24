from libcpp cimport bool as cppbool
from libcpp.vector cimport vector as cppvector
from gmpy2 cimport mpz 

cdef extern from "gmp.h":
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
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)

ctypedef mp_limb_t ulong
ctypedef mp_limb_signed_t slong

cdef extern from "flint/fmpz.h":
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]
    void fmpz_get_mpz(mpz_t x, const fmpz_t f) noexcept
    void fmpz_init(fmpz_t f)
    void fmpz_set_ui(fmpz_t f, ulong g)
    void fmpz_set_si(fmpz_t f, slong g)
    void fmpz_clear(fmpz_t f)
    int fmpz_print(const fmpz_t x)
cdef extern from "pplite/pplite.hh" namespace "pplite":
# cdef extern from "pplite/FLINT_Integer.hh" namespace "pplite":
    cdef cppclass FLINT_Integer
    cdef cppclass FLINT_Integer:
        FLINT_Integer()
        const fmpz_t& impl() 
        FLINT_Integer(const fmpz_t z)
        FLINT_Integer(signed int si)
        FLINT_Integer(const mpz_t z)
# cdef extern from "pplite/GMP_Integer.hh" namespace "pplite":
#     cdef cppclass GMP_Integer
#     cdef cppclass GMP_Integer:
#         GMP_Integer()
#         const mpz_t& impl() 
#         # GMP_Integer(const fmpz_t z)
#         GMP_Integer(signed int si)
#         GMP_Integer(const mpz_t z)



# cdef extern from "pplite/globals.hh" namespace "pplite":
    ctypedef size_t dim_type

# cdef extern from "pplite/Var.hh" namespace "pplite":
    cdef cppclass Var
    cdef cppclass Var:
#       Var()
        Var(dim_type i)
        dim_type id()
        dim_type space_dim()
        # void m_swap(Var& v)
    ctypedef struct Vars_Set
    cppbool less_than(Var v, Var w)
    void swap(Var& v, Var& w)
    
# cdef extern from "pplite/Linear_Expr.hh" namespace "pplite":
    ctypedef cppvector[FLINT_Integer] Impl
    cdef cppclass Linear_Expr
    cdef cppclass Linear_Expr:
        Linear_Expr()
        Linear_Expr(Var v)
        Linear_Expr(dim_type dim)
        Linear_Expr(Linear_Expr &e)
        Linear_Expr(const Linear_Expr &e, dim_type dim)
        Linear_Expr operator+(Linear_Expr &e)
        # Linear_Expr operator+(Var v, Linear_Expr e)
        Linear_Expr operator-(Linear_Expr &e) # which is correct?
        Linear_Expr operator*(FLINT_Integer &c)
        # Con operator >(Linear_Expr &e)
        # Con operator <(Linear_Expr &e)
        # Con operator <=(Linear_Expr &e)
        # Con operator >=(Linear_Expr &e)
        # Con operator ==(Linear_Expr &e)       
        dim_type id()  #methods? #linear_expr.hh lines 39-112
        dim_type space_dim()
        FLINT_Integer get(dim_type dim)
        FLINT_Integer get(Var v)
        void set(dim_type dim, FLINT_Integer n)
        void set(Var v, FLINT_Integer n)
        void set_space_dim(dim_type dim)
        void swap_space_dims(dim_type i, dim_type j)
        void shift_space_dims(dim_type start, dim_type n)
        void shift_space_dims(Var start, dim_type n)
        # void remove_space_dims(Iter first, Iter last)
        Impl impl()
        cppbool is_equal_to(Linear_Expr& y) 
        cppbool is_zero()
    # Linear_Expr operator+(Linear_Expr e1, Linear_Expr e2)
    Linear_Expr operator+(Linear_Expr e, Var v) 
    Linear_Expr operator+(Var v, Linear_Expr v)
    Linear_Expr operator+(Var v, Var w)
    # Linear_Expr operator-(Linear_Expr e1, Linear_Expr e2)
    Linear_Expr operator-(Linear_Expr &e, Var &v)
    Linear_Expr operator-(Var &v, Linear_Expr &v)
    Linear_Expr operator-(Var v, Var w)
    void neg_assign(Linear_Expr& e)


# cdef extern from "pplite/Affine_Expr.hh" namespace "pplite":
    cdef cppclass Affine_Expr
    cdef cppclass Affine_Expr:
        Linear_Expr expr
        FLINT_Integer inhomo
        Affine_Expr()
        Affine_Expr(FLINT_Integer i)
        Affine_Expr(Linear_Expr e, FLINT_Integer i)
        Affine_Expr(Affine_Expr &e)
        dim_type space_dim()
        void set_space_dim(dim_type dim)
        cppbool is_zero()
        void m_swap(Affine_Expr& y)
        void normalize()
        void sign_normalize()
        Affine_Expr operator+(Affine_Expr &a)
    #    Affine_Expr operator+(Linear_Expr &e, FLINT_Integer &n)
        Affine_Expr operator-(Affine_Expr &a)
        Affine_Expr operator*(FLINT_Integer &c)
    # Affine_Expr operator+(Linear_Expr &e, Var &v)
    # Affine_Expr operator+(Linear_Expr &e, FLINT_Integer &n)
    # Affine_Expr operator+(Affine_Expr &e, FLINT_Integer &n)
    # Affine_Expr operator+(Affine_Expr &a1, Linear_Expr &e2)
    # Affine_Expr operator+(Affine_Expr &a, Linear_Expr &e2)
    Affine_Expr operator+(Affine_Expr a, Var v)
    Affine_Expr operator-(Var v, Affine_Expr a)
    Affine_Expr operator-(Linear_Expr e1, Affine_Expr a1)
    void neg_assign(Affine_Expr& a)
    # Affine_Expr& operator+=(Affine_Expr& a1, Var v) #+= operator not yet supported. FML.
# cdef extern from "pplite/Con.hh":
#     ctypedef enum ConType "Con::Type":
#         EQUALITY, NONSTRICT_INEQUALITY, STRICT_INEQUALITY        
    cdef cppclass Con
    cdef cppclass Con:
        enum ConType "Con::Type":
            EQUALITY
            NONSTRICT_INEQUALITY
            STRICT_INEQUALITY
        struct Impl:
            Linear_Expr expr
            FLINT_Integer inhomo
            ConType type
            Impl()
            Impl(Linear_Expr e, FLINT_Integer i, ConType t)
        Con()
        Con(const Con &c)
        Con(Linear_Expr expr, FLINT_Integer inhomo, ConType type)
        Con(Affine_Expr ae, ConType type)
        dim_type space_dim()    
        void set_space_dim(dim_type dim)
        # void permute_space_dims_cycle(const Dims& cycle, dim_type d)
        Impl& impl()
        ConType  type()
        cppbool is_equality()
        cppbool is_inequality()
        cppbool is_nonstrict_inequality()
        cppbool is_strict_inequality()
        Linear_Expr linear_expr()
        FLINT_Integer coeff(Var v)
        FLINT_Integer inhomo_term()
        Con zero_dim_false()
        

        # con.hh lines 215-489. 
    Con operator=(Con &c)
    Con operator<(Linear_Expr e1, const Linear_Expr& e2)
    Con operator<(Var v1, Var v2)
    Con operator<(Linear_Expr e, const FLINT_Integer& n)
    Con operator<(FLINT_Integer n, Linear_Expr e)
    Con operator>(Linear_Expr e1, const Linear_Expr& e2)
    Con operator>(Var v1, Var v2)
    Con operator>(Linear_Expr e, FLINT_Integer n)
    Con operator>(FLINT_Integer n, Linear_Expr e)
    Con operator==(Linear_Expr e1, const Linear_Expr& e2)
    Con operator==(Var v1, Var v2)          
    Con operator==(Linear_Expr e, FLINT_Integer n)
    Con operator==(FLINT_Integer n, Linear_Expr e)
    Con operator<=(Linear_Expr e1, const Linear_Expr& e2)
    Con operator<=(Var v1, Var v2)
    Con operator<=(Linear_Expr e, FLINT_Integer n)
    Con operator<=(FLINT_Integer n, Linear_Expr e)
    Con operator>=(Linear_Expr e1, const Linear_Expr& e2)
    Con operator>=(Var v1, Var v2)
    Con operator>=(Linear_Expr e, FLINT_Integer n)
    Con operator>=(FLINT_Integer n, Linear_Expr e)
    Con operator<(Affine_Expr e1, const Affine_Expr& e2)
    Con operator<(Affine_Expr e, const FLINT_Integer& n)
    Con operator<(Affine_Expr e, Var v)
    Con operator<(const FLINT_Integer& n, Affine_Expr e)
    Con operator<(Var v, Affine_Expr e)
    Con operator>(Affine_Expr e1, const Affine_Expr& e2)
    Con operator>(Affine_Expr e, const FLINT_Integer& n)
    Con operator>(Affine_Expr e, Var v)
    Con operator>(const FLINT_Integer& n, Affine_Expr e)
    Con operator>(Var v, Affine_Expr e)
    Con operator==(Affine_Expr e1, const Affine_Expr& e2)
    Con operator==(Affine_Expr e, const FLINT_Integer& n)
    Con operator==(Affine_Expr e, Var v)
    Con operator==(const FLINT_Integer& n, Affine_Expr e)
    Con operator==(Var v, Affine_Expr e)
    Con operator<=(Affine_Expr e1, const Affine_Expr& e2)
    Con operator<=(Affine_Expr e, const FLINT_Integer& n)
    Con operator<=(Affine_Expr e, Var v)    
    Con operator<=(const FLINT_Integer& n, Affine_Expr e)
    Con operator<=(Var v, Affine_Expr e)
    Con operator>=(Affine_Expr e1, const Affine_Expr& e2)
    Con operator>=(Affine_Expr e, const FLINT_Integer& n)
    Con operator>=(Affine_Expr e, Var v)
    Con operator>=(const FLINT_Integer& n, Affine_Expr e)
    Con operator>=(Var v, Affine_Expr e)

    # Con operator<(Affine_Expr e1, const Affine_Expr& e2)
    # Con operator>(Affine_Expr e1, const Affine_Expr& e2)
    # Con operator<=(Affine_Expr e1, const Affine_Expr& e2)
    # Con operator>=(Affine_Expr e1, const Affine_Expr& e2)
    #     # Con operator!=(Affine_Expr e1, const Affine_Expr& e2)
    # Con operator==(Affine_Expr e1, const Affine_Expr& e2)
    #     # Con operator
    # Con operator<(Linear_Expr e1, Linear_Expr e2)
    # Con operator<(Var e1, Var e2)
    # Con operator<(Var v1, Linear_Expr e2)
    # Con operator<(Var e1, Var e2)
    # Con operator>(Linear_Expr e1, Linear_Expr e2)
    # Con operator>(Var e1, Var e2)
    # Con operator<=(Linear_Expr e1, Linear_Expr e2)
    # Con operator<=(Var e1, Var e2)
    # Con operator>=(Linear_Expr e1, Linear_Expr e2)
    # Con operator>=(Var e1, Var e2)
    # Con operator==(Linear_Expr e1, Linear_Expr e2)
    # Con operator==(Var e1, Var e2)
        # Con operator!=(Linear_Expr e1, Linear_Expr e2)
        # Con operator!=(Var e1, Var e2)  

