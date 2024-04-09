# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check
# from flint cimport fmpz
from libcpp.vector cimport vector as cppvector
import_gmpy2()

cdef class Flint_Int(object):
    r"""
    Wrapper for FLINT_Integer
    
    Stuff
    
    """
    def __cinit__(self, int z):
        cdef fmpz_t fmpz_int
        fmpz_init(fmpz_int)
        fmpz_set_ui(fmpz_int, z)
        self.thisptr = new FLINT_Integer(fmpz_int)
        fmpz_clear(fmpz_int)
    def __dealloc__(self):
        del self.thisptr
    def __repr__(self):
        self.thisptr.print()

# cdef FLINT_Integer_to_python(FLINT_Integer integer):
#     cdef mpz_t new_integer
#     fmpz_get_mpz(new_integer, integer.impl())
#     py_int = GMPy_MPZ_From_mpz(new_integer)
#     mpz_clear(new_integer)
#     return py_int



@cython.freelist(128)
cdef class Variable(object):
    r"""
    Wrapper for PPL's ``Variable`` class.

    A dimension of the vector space.    
    
    OTHER STUFF HERE

    INPUT:

    - ``i`` -- integer. The index of the axis.

    OUTPUT:

    A :class:`Variable`

    Examples:

    >>> from ppl import Variable
    >>> x = Variable(123)
    >>> x.id()
    123
    >>> x
    x123

    Note that the "meaning" of an object of the class Variable is completely
    specified by the integer index provided to its constructor: be careful not
    to be mislead by C++ language variable names. For instance, in the following
    example the linear expressions ``e1`` and ``e2`` are equivalent, since the
    two variables ``x`` and ``z`` denote the same Cartesian axis:

    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> z = Variable(0)
    >>> e1 = x + y; e1
    x0+x1
    >>> e2 = y + z; e2
    x0+x1
    >>> e1 - e2
    0
    """
    def __cinit__(self, dim_type i):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        Tests:

        >>> from pplite import Variable
        >>> Variable(123)   # indirect doctest
        x123
        """
        self.thisptr = new Var(i)

        
    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr
        
    def __hash__(self):
        r"""
        Tests:

        >>> import pplite
        >>> hash(ppl.Variable(12))
        Traceback (most recent call last):
        ...
        TypeError: Variable unhashable
        """
        raise TypeError('Variable unhashable')

    def id(self):
        """
        Return the index of the Cartesian axis associated to the variable.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(123)
        >>> x.id()
        123
        """
        return self.thisptr.id()
        
    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer. The returned value is ``self.id()+1``.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> x.space_dimension()
        1
        """
        return self.thisptr.space_dim()
    
    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from pplite import Variable
        >>> x = Variable(0)
        >>> x.__repr__()
        'x0'
        """
        return 'x{0}'.format(self.id())
    
    def __add__(self, other):
        r"""
        Return the sum ``self`` + ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` + ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x + 15
        x0+15
        >>> 15 + y
        x1+15

        >>> from gmpy2 import mpz
        >>> x + mpz(3)
        x0+3
        >>> mpz(-5) + y
        x1-5

        >>> x + 1.5
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        >>> 1.5 + x
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        """
        return Linear_Expression(self) + Linear_Expression(other)

    def __radd__(self, other):
        return Linear_Expression(self) + Linear_Expression(other)

    def __sub__(self, other):
        r"""
        Return the difference ``self`` - ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` - ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x - 15
        x0-15
        >>> 15 - y
        -x1+15
        """
        return Linear_Expression(self) - Linear_Expression(other)

    def __rsub__(self, other):
        return Linear_Expression(other) - Linear_Expression(self)

    def __mul__(self, other):
        r"""
        Return the product ``self`` * ``other``.

        INPUT:

        - ``self``, ``other`` -- One must be an integer, the other a
          :class:`Variable`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` * ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x * 15
        15*x0
        >>> 15 * y
        15*x1

        >>> 1.5 * x
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        >>> x * 1.5
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        """
        if isinstance(self, Variable):
            return Linear_Expression(self) * other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            return Linear_Expression(other) * self

    def __rmul__(self, other):
        return Linear_Expression(self) * other

    def __pos__(self):
        r"""
        Return ``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``+self``

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); x
        x0
        >>> +x
        x0
        """
        return Linear_Expression(self)

    def __neg__(self):
        r"""
        Return -``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``-self``

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); x
        x0
        >>> -x
        -x0
        """
        return Linear_Expression(self)*(-1)
    

####################################################
### Linear_Expression ##############################
####################################################
cdef class Linear_Expression(object):
    r"""
    Wrapper for PPL's ``PPL_Linear_Expression`` class.
    """
    def __init__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        Tests:

        >>> from ppl import Linear_Expression
        >>> Linear_Expression(10)   # indirect doctest
        10
        """
        cdef size_t i
        if len(args) == 2:
            a = args[0]
            b = args[1]
            # self.thisptr = new Linear_Expr()
            # if isinstance(a, dict):
            #     if a:
            #         self.thisptr.set_space_dimension(1 + max(a))
            #         for i, coeff in a.items():
            #             self.thisptr.set_coefficient(Var(i), PPL_Coefficient_from_pyobject(coeff))
            # else:
            #     self.thisptr.set_space_dimension(len(a))
            #     for i, coeff in enumerate(a):
            #         self.thisptr.set_coefficient(PPL_Variable(i), PPL_Coefficient_from_pyobject(coeff))
            # self.thisptr.set_inhomogeneous_term(PPL_Coefficient_from_pyobject(b))
            # return
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, Variable):
                v = <Variable>arg
                self.thisptr = new Linear_Expr(v.thisptr[0])
                return
            if isinstance(arg, Linear_Expression):
                e = <Linear_Expression>arg
                self.thisptr = new Linear_Expr(e.thisptr[0])
                return
            self.thisptr = new Linear_Expr()
        elif len(args) == 0:
            self.thisptr = new Linear_Expr()
            return
        else:
            raise ValueError("Cannot initialize with more than 2 arguments.")

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Linear_Expression(10))
        Traceback (most recent call last):
        ...
        TypeError: Linear_Expression unhashable
        """
        raise TypeError('Linear_Expression unhashable')

    def space_dimension(self):
        """
        Return the dimension of the vector space necessary for the
        linear expression.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x+y+1).space_dimension()
        2
        >>> (x+y).space_dimension()
        2
        >>> (y+1).space_dimension()
        2
        >>> (x+1).space_dimension()
        1
        >>> (y+1-y).space_dimension()
        2
        """
        return self.thisptr.space_dim()

#     def get_coeffients_list(self):
#         # Problem, how to get python to access the data contained in the class construction. These methods aren't written explicltily or done for me in friendly way.
        
#         # Goal. Convert Impl being a vecotr of integers to a list of integers in python or to a python object or be able to access particular elements of the vector Impl as python objects.
        
#         # Impl is an iterating object so I should be able to iterate over it. 
        
#         # Last time: Impl, standard vector. This is the place in the class decleartions where the data is stored. These are private members and can be access via the public method of the impl() method. Assuming that I have read this correctly.
        
#         # Impl is a vector (in the standard libary) of the type Integer. 
#         # The type Integer might be defined defined in Integer_fwd.hh. I am unsure where this type is defined. 
#         # The integer type is used to alais the integers used for computaions - either FLINT or GMP. 
#         #
        
#         # So where does Integer come from?
#         Other questions include:
#         How to take in new data because in the constructor,   explicit Linear_Expr(dim_type dim) : row(dim)

        
    def coefficient(self, v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An Integer. 

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> e = 3*x+1
        >>> e.coefficient(x)
        mpz(3)
        >>> e.coefficient(Variable(13))
        mpz(0)
        """
        cdef Variable vv
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        #return FLINT_Integer_to_python(self.thisptr.impl()[vv.id()])
#     def inhomogeneous_term(self):
#         """
#         Return the inhomogeneous term of the linear expression.

#         OUTPUT:

#         Integer.

#         Examples:

#         >>> from ppl import Linear_Expression
#         >>> Linear_Expression(10).inhomogeneous_term()
#         mpz(10)
#         """
#         return GMPy_MPZ_From_mpz(self.thisptr.inhomogeneous_term().get_mpz_t())
    
#     def __repr__(self):
#         r"""
#         Return a string representation of the linear expression.

#         OUTPUT:

#         A string.

#         Examples:

#         >>> from ppl import Linear_Expression, Variable
#         >>> x = Variable(0)
#         >>> y = Variable(1)
#         >>> x+1
#         x0+1
#         >>> x+1-x
#         1
#         >>> 2*x
#         2*x0
#         >>> x-x-1
#         -1
#         >>> x-x
#         0
#         """
#         s = ''
#         first = True
#         for i in range(self.space_dimension()):
#             x = Variable(i)
#             coeff = self.coefficient(x)
#             if coeff == 0:
#                 continue
#             if first and coeff == 1:
#                 s += '%r' % x
#                 first = False
#             elif first and coeff == -1:
#                 s += '-%r' % x
#                 first = False
#             elif first and coeff != 1:
#                 s += '%d*%r' % (coeff, x)
#                 first = False
#             elif coeff == 1:
#                 s += '+%r' % x
#             elif coeff == -1:
#                 s += '-%r' % x
#             else:
#                 s += '%+d*%r' % (coeff, x)
#         inhomog = self.inhomogeneous_term()
#         if inhomog != 0:
#             if first:
#                 s += '%d' % inhomog
#                 first = False
#             else:
#                 s += '%+d' % inhomog
#         if first:
#             s = '0'
#         return s



def test_current_obj():
    # cdef fmpz_t x
    # fmpz_init(x)
    # fmpz_set_ui(x, 7)
    # w = FLINT_Integer(x)
    # fmpz_clear(x)
    # z = FLINT_Integer_to_python(w)
    z = Flint_Int(7)
    print(z)
    # x = Variable(10)
    # e = Linear_Expression(x)
    # x_2 = Variable(2)
    # print(e.coefficient(x))
    # print(e.coefficient(x_2))