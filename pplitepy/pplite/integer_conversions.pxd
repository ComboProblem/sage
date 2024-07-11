from __future__ import absolute_import

from .pplite_decl cimport *
####################################################

cdef FLINT_Integer_to_Python(FLINT_Integer& integer)

cdef FLINT_Integer Python_int_to_FLINT_Integer(integer)