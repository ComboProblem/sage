# distutils: language = c++
# distutils: libraries = gmp gmpxx pplite m flint

cimport cython

from gmpy2 cimport import_gmpy2, mpz, mpz_t, GMPy_MPZ_From_mpz, MPZ_Check, mpq, MPQ_Check, GMPy_MPQ_From_mpz
from libcpp.vector cimport vector as cppvector

# from numbers import Rational

#  about inclusion of types. concrete type to handle is faraciton. 

import_gmpy2()

cdef FLINT_Integer_to_Python(FLINT_Integer& integer):
    """ Converts FLINT_Integer to python object.

    INPUT: :class:`FLINT_Interger` (c++)
    OUTPUT: Python object
    """
    cdef mpz_t new_int
    mpz_init(new_int)
    fmpz_get_mpz(new_int, integer.impl())
    y = GMPy_MPZ_From_mpz(new_int)
    mpz_clear(new_int)
    return y

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer):
    cdef fmpz_t x
    cdef fmpz y
    if isinstance(integer, (int, str)):
        fmpz_init(x)
        fmpz_set_si(x, integer)
    return FLINT_Integer(x)
    if MPZ_Check(integer): # is this okay?
        y = <fmpz> integer
        return FLINT_Integer(y)
    raise ValueError("Integer Conversion Failed")


cdef FLINT_Rational_to_Python(FLINT_Rational& rational):
    """Converts the Flint_Rational class to a python object."""
    cdef mpz_t a
    cdef mpz_t b
    mpz_init(a)
    mpz_init(b)
    fmpq_get_mpz_frac(a , b, rational.impl())
    frac = GMPy_MPQ_From_mpz(a, b)
    mpz_clear(a)
    mpz_clear(b)
    return frac

cdef FLINT_Rational Python_float_to_FLINT_Rational(rational):
    """ Converts python float or fraction """
    cdef FLINT_Integer num
    cdef FLINT_Integer den
    try:
        numerator, denominator = rational.as_integer_ratio()
    except ValueError:
        raise ValueError("Rational Conversion Failed.")
    num = Python_int_to_FLINT_Integer(numerator)
    dem = Python_int_to_FLINT_Integer(denominator)
    return FLINT_Rational(num, dem)




def FLINT_Integer_Conversion_Check(possible_integer):
    """
    Checks a python object is convertible to a FLINT_Integer.

    Input: Object

    Output: Bool
    """
    if isinstance(possible_integer, (int, str)):
        return True
    if MPZ_Check(possible_integer):
        return True
    return False