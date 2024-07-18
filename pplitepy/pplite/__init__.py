r"""
Cython wrapper for the Parma Polyhedra Lite Library (PPLite)

"""

__version__ = "0.0.26"

from .linear_algebra import (
        Variable, Linear_Expression, Affine_Expression
        )

from .constraint import (
        Constraint
        )

from .generators import (
        PPliteGenerator
        )
