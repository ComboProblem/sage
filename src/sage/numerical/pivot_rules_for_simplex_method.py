r"""
Pivot Rules for Simplex Method

This module is an addendum to the Interactive Simplex Method
to draw connection between active research on simplex method pivot rules
and an educational settings.

A pivot rule is a function or algorithm on an linear program (LP) `\min\{cx: Ax<=b\}`
- (a) `v` is adjacent to `u` in the graph of $P$ and
- (b) `cu \leq cv`.

The module offers a number of predefined pivot rules and classes which serve as plug and play for defining/experimenting with pivot rules.

Authors:

-Acadia Larsen (March 1, 2024): Initial Version.

EXAMPLES:

This module is designed for use with :mod:`sage.interactive_simplex_method`::
    sage: A = ([1, 1], [3, 1], [-1, -1])
    sage: b = (1000, 1500, -400)
    sage: c = (10, 5)
    sage: P = InteractiveLPProblemStandardForm(A, b, c)
    sage: P.run_simplex_method('steepest_edge')
    \begin{equation*}
    ...
    \end{equation*}
    The initial dictionary is infeasible, solving auxiliary problem.
    ...
    Entering: $x_{0}$. Leaving: $x_{5}$.
    ...
    Entering: $x_{1}$. Leaving: $x_{0}$.
    ...
    Back to the original problem.
    ...
    Entering: $x_{5}$. Leaving: $x_{4}$.
    ...
    Entering: $x_{2}$. Leaving: $x_{3}$.
    ...
    The optimal value: $6250$. An optimal solution: $\left(250,\,750\right)$.

It allows one to explore different pivot rules such as the Normalized Weight pivot rule ([BDLS2022]_.)::

    sage: A = ([1, 1], [3, 1], [-1, -1])
    sage: b = (1000, 1500, -400)
    sage: c = (10, 5)
    sage: P = InteractiveLPProblemStandardForm(A, b, c)
    sage: def eta(v):
    ....:     return v.norm()
    sage: D = P.auxiliary_problem().initial_dictionary() # Phase I
    sage: D.enter('x0')
    sage: D.leave('x5')
    sage: D.update()
    sage: D.run_simplex_method('NW_rule', eta, [1,1,4])
    \begin{equation*}
    ...
    \end{equation*}
    Entering: $x_{2}$. Leaving: $x_{0}$.
    \begin{equation*}
    ...
    \end{equation*}
    sage: d = P.feasible_dictionary(D) # Phase II
    sage: d.run_simplex_method('NW_rule', eta, [1,1])
    \begin{equation*}
    ...
    \end{equation*}
    Entering: $x_{5}$. Leaving: $x_{3}$.
    \begin{equation*}
    ...
    \end{equation*}
    Entering: $x_{1}$. Leaving: $x_{4}$.
    \begin{equation*}
    ...
    \end{equation*}

Pivot Rules can be used to manipulate dictionaries defined by an :class:`~sage.src.numerical.interactive_simplex_method.LPDictionary` on a step by step basis::

    sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
    sage: A = ([1, 1], [3, 1])
    sage: b = (1000, 1500)
    sage: c = (10, 5)
    sage: P = InteractiveLPProblemStandardForm(A, b, c)
    sage: D = P.initial_dictionary()
    sage: pivot = SimplexMethodPivot("dantzig")
    sage: pivot(D)
    sage: D.entering()
    x1
    sage: D.leaving()
    x4

Define custom pivot rules::

    sage: from sage.numerical.pivot_rules_for_simplex_method import AbstractSimplexMethodPivotRule
    sage: class MyPivotRule(AbstractSimplexMethodPivotRule):
    ....:     def dictionary_pivot(dict):
    ....:         dict.enter('x1')
    ....:         dict.leave('x4')
    sage: my_pivot = SimplexMethodPivot(MyPivotRule)
    sage: A = ([1, 1], [3, 1])
    sage: b = (1000, 1500)
    sage: c = (10, 5)
    sage: P = InteractiveLPProblemStandardForm(A, b, c)
    sage: D = P.initial_dictionary()
    sage: my_pivot(D)
    sage: D.entering()
    x1
    sage: D.leaving()
    x4

Classes and functions
---------------------
"""
# ****************************************************************************
#       Copyright (C) 2024 Acadia Larsen acadia.larsen@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import deepcopy
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.all import SageObject
from sage.misc.abstract_method import abstract_method
from sage.modules.free_module_element import free_module_element as vector
from sage.rings.infinity import InfinityRing
from numpy import argmax
#from sage.numerical.interactive_simplex_method import LPAbstractDictionary

class AbstractSimplexMethodPivotRule(UniqueRepresentation):
    r"""
    Abstract class for memoryless pivot rules.

    Individual pivot rules are not in Sage's namespace.
    Using a pivot rule should be instantiated via :class:`~sage.numerical.pivot_rules_for_simplex_method.SimplexMethodPivot`.
    """
    @classmethod
    def has_dual_method(self):
        r""" For use in SimplexMethodPivot to identify if a dual pivot rule has been written.

        Overwrite this method after implementing a dual pivot method. See ``BlandsRule`` for example.

        OUTPUT: bool
        """
        return False
    @abstract_method
    def dictionary_pivot(dictionary, *args):
        r"""
        Abstract method for computing a pivot for the primal simplex method based on an :class:`sage.numerical.interactive_simplex_method.LPDictionary`.

        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` for the primal simplex method.

        INPUT:

        - ``dictionary`` -- an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`

        - ``*args`` -- required inputs for specific method

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import BlandsRule
            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.entering()

            sage: D.leaving()

            sage: BlandsRule.dictionary_pivot(D)
            sage: D.entering()
            x1
            sage: D.leaving()
            x4

        Handling unboundedness, feasibility, and optimality should be left to the solver.::

            sage: A = ([1, 0],)
            sage: b = (1,)
            sage: c = (0, 1)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.run_simplex_method()
            \begin{equation*}
            ...
            \end{equation*}
            The problem is unbounded in $x_{2}$ direction.

        The leaving variable should be ``None`` if a dictionary represents unboundedness in a problem::

            sage: D = P.initial_dictionary()
            sage: BlandsRule.dictionary_pivot(D)
            sage: D.entering()
            x2
            sage: D.leaving() is None
            True

        The :meth:`__call__` method of :class:`SimplexMethodPivot` is responsible checking for type input::

            sage: LP_tablaux = Matrix([[0, 10, 5], [1, 1, 1000], [3, 1, 1500]])
            sage: BlandsRule.dictionary_pivot(LP_tablaux)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matrix.matrix_integer_dense.Matrix_integer_dense' object has no attribute 'entering'

            """
    @abstract_method
    def dual_dictionary_pivot(dictionary, *args):
        r"""
        Abstract method for computing a pivot for the dual simplex method based on an :class:`sage.numerical.interactive_simplex_method.LPDictionary`.

        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` for the dual simplex method.

        One does not need to define a :meth:`dual_dictionary_pivot` unless intending to run a dual simplex algorithm.

        INPUT:

        - ``dictionary`` -- an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`.

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import BlandsRule
            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary(2, 3, 5)
            sage: D.entering()

            sage: D.leaving()

            sage: BlandsRule.dual_dictionary_pivot(D)
            sage: D.entering()
            x1
            sage: D.leaving()
            x3

        Some pivot rules do not have dual methods implemented::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SteepestEdge
            sage: SteepestEdge.has_dual_method()
            False
        """

class SimplexMethodPivot(SageObject, UniqueRepresentation):
    r"""
    A class to perform pivots in the simplex method of solving Linear Programs.

    INPUT:

    - ``name`` -- a string or a ``SimplexMethodPivotRule`` object.

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: blands_rule = SimplexMethodPivot()
        sage: D = P.initial_dictionary()
        sage: D.entering()

        sage: D.leaving()

        sage: blands_rule(D)
        sage: D.entering()
        x1
        sage: D.leaving()
        x4
    """
    @staticmethod
    def __classcall__(cls, name=None):
        r"""
        Normalize init agrugments.

        INPUT:

        - ``name`` -- a string or a ``SimplexMethodPivotRule`` object

        TESTS::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: from sage.numerical.pivot_rules_for_simplex_method import BlandsRule
            sage: SimplexMethodPivot("blands_rule") is SimplexMethodPivot()
            True
            sage: SimplexMethodPivot("blands_rule") is SimplexMethodPivot(BlandsRule)
            True
        """
        if name == "blands_rule" or name is None:
            return super().__classcall__(cls, pivot_rule=BlandsRule)
        if name == "steepest_edge":
            return super().__classcall__(cls, pivot_rule=SteepestEdge)
        if name == "dantzig":
            return super().__classcall__(cls, pivot_rule=DantzigsRule)
        if name == "NW_rule":
            return super().__classcall__(cls, pivot_rule=NWRule)
        if issubclass(name, AbstractSimplexMethodPivotRule):
            return super().__classcall__(cls, pivot_rule=name)
        else:
            raise TypeError("Input is required to be a string with the name of a predefined pivot rule or an instance of ``AbstractSimplexMethodPivotRule``")

    def __init__(self, pivot_rule):
        r"""
        Initialize a pivot rule. This is never directly called.

        INPUT:

        - ``pivot_rule`` -- a ``SimplexMethodPivotRule``

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot

        Create Blands Rule::

            sage: blands_rule = SimplexMethodPivot("blands_rule")

        No input of a named pivot rule uses Blands Rule instead::

            sage: generic_pivot_rule = SimplexMethodPivot()
            sage: blands_rule is generic_pivot_rule
            True

        Create Steepest Edge Pivot Rule::

            sage: steepest_edge = SimplexMethodPivot("steepest_edge")

        Create Dantzigs Pivot Rule::

            sage: dantzigs_rule =  SimplexMethodPivot("dantzig")

        Create Normalized Weight Pivot Rule::

            sage: NW_rule =  SimplexMethodPivot("NW_rule")
        """
        self._pivot_rule = pivot_rule
        super().__init__()

    def __call__(self, dictionary, dual_mode=False, *args):
        r"""
        Makes pivot rules have a function like notation and normalizes inputs.

        INPUT:

        - ``dictionary`` -- instance of :LPDictionary:

        - ``args`` --  arguments for pivot rule

        EXAMPLES::

            A pivot rule applied to an :class:`sage.numerical.interactive_simplex_method.LPDictionary` sets the entering
            and leaving variables::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: blands_rule = SimplexMethodPivot()
            sage: D.entering()

            sage: D.leaving()

            sage: blands_rule(D)

            sage: D.entering()
            x1
            sage: D.leaving()
            x4

        The call method also checks types::

            sage: LP_tablaux = Matrix([[0, 10, 5], [1, 1, 1000], [3, 1, 1500]])
            sage: blands_rule = SimplexMethodPivot()
            sage: blands_rule(LP_tablaux)
            Traceback (most recent call last):
            ...
            TypeError: Input is required to be an instance of `LPAbstractDictionary`
        """
        if dual_mode:
            if self._pivot_rule.has_dual_method():
                pivot = self._pivot_rule.dual_dictionary_pivot
            else:
                raise NotImplementedError("The selected pivot rule does not have an implementation that is compatible with the dual simplex method.")
        else:
            pivot = self._pivot_rule.dictionary_pivot
        if isinstance(dictionary, LPAbstractDictionary):
            if args:
                pivot(dictionary, *args)
            else:
                pivot(dictionary)
        else:
            raise TypeError("Input is required to be an instance of `LPAbstractDictionary`")

class BlandsRule(AbstractSimplexMethodPivotRule):
    r"""
    Defines Blands Rule as a ``AbstractSimplexMethodPivotRule``.
    Ref. [TZ1993]_

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: D.entering()

        sage: D.leaving()

        sage: pivot = SimplexMethodPivot("blands_rule")
        sage: pivot(D)
        sage: D.entering()
        x1
        sage: D.leaving()
        x4
    """
    def has_dual_method():
        return True

    def dictionary_pivot(current_dictionary):
        r"""
        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` as defined by Blands rule.

        INPUT:

        - ``current_dictionary`` -- an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: pivot = SimplexMethodPivot("blands_rule") #indirect test
            sage: pivot(D)
            sage: D.entering()
            x1
            sage: D.leaving()
            x4
        """
        if current_dictionary.entering() is None:
            entering_variable = min(current_dictionary.possible_entering())
            current_dictionary.enter(entering_variable)
        if current_dictionary.leaving() is None:
            possible = current_dictionary.possible_leaving()
            if possible:
                leaving_variable = min(possible)
                current_dictionary.leave(leaving_variable)

    def dual_dictionary_pivot(current_dictionary):
        r"""
        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` as defined by Blands rule for dual simplex method.

        INPUT:

        - ``current_dictionary`` --  an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary(2, 3, 5) # dual feasible dictionary
            sage: pivot = SimplexMethodPivot()
            sage: pivot(D, True)
            sage: D.leaving()
            x3
            sage: D.entering()
            x1
        """
        if current_dictionary.leaving() is None:
            entering_variable = min(current_dictionary.possible_leaving())
            current_dictionary.leave(entering_variable)
        if current_dictionary.entering() is None:
            possible = current_dictionary.possible_entering()
            if possible:
                leaving_variable = min(possible)
                current_dictionary.enter(leaving_variable)

class SteepestEdge(AbstractSimplexMethodPivotRule):
    r"""
    Defines the steepest edge pivot rule as a :class:`AbstractSimplexMethodPivotRule`.
    The entering variable is selected by finding the steepest edge.
    The leaving variable is selected by finding the minimal ratio of the possible leaving variables.
    This method is generalized to selection based on the `p` norm with `p \geq 1`.
    For more information see [GR1977]_ and [KMT2019]_.

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: D.entering()

        sage: D.leaving()

        sage: pivot = SimplexMethodPivot("steepest_edge")
        sage: pivot(D)
        sage: D.entering()
        x1
        sage: D.leaving()
        x4
    """

    def dictionary_pivot(current_dictionary, p=2):
        r"""
        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` as defined by steepest edge rule.

        Input:

        - ``current_dictionary`` -- an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`

        - ``p`` -- optional real number with `p\\geq 1`

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: pivot = SimplexMethodPivot("steepest_edge") # indirect test
            sage: pivot(D, False, 1) # let's use the 1 norm.
            sage: D.entering()
            x1
            sage: D.leaving()
            x4

            Values of p are checked::

            sage: D.update()
            sage: pivot(D, False, 0)
            Traceback (most recent call last):
            ...
            ValueError: p must be at least 1
        """
        if p < 1:
            raise ValueError("p must be at least 1")
        if current_dictionary.entering() is None:
            thetas = []
            for v in current_dictionary.basic_variables():
                # We want to avoid p-root complications, so work with the pth power.
                thetas.append((1+sum([i**p for i in current_dictionary.row_coefficients(v)])))
            best_value = -1*current_dictionary.objective_coefficients()[0]/thetas[0]
            index_to_select = 0
            # This loop selects the index that corresponds to the steepest edge.
            for i in range(1,len(current_dictionary.nonbasic_variables())):
                decent_direction = -1*current_dictionary.objective_coefficients()[i]/thetas[i-1]
                if decent_direction < best_value:
                    best_value = decent_direction
                    index_to_select = i
            entering_variable = current_dictionary.nonbasic_variables()[index_to_select]
            current_dictionary.enter(entering_variable)
        if current_dictionary.leaving() is None:
            possible = current_dictionary.possible_leaving()
            if possible:
                # We select leaving variables based on ratios.
                leaving_variable = min(current_dictionary.ratios())[1]
                current_dictionary.leave(leaving_variable)

class DantzigsRule(AbstractSimplexMethodPivotRule):
    r"""
    Selects the entering variable to be the variable with the most negative reduced cost.

    Originally proposed by Dantzig in the formulation of linear programming. Also know as greatest improvement rule.

    See [TZ1993]_.

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: D.entering()

        sage: D.leaving()

        sage: pivot = SimplexMethodPivot("dantzig")
        sage: pivot(D)
        sage: D.entering()
        x1
        sage: D.leaving()
        x4
    """
    def dictionary_pivot(current_dictionary):
        r"""
        Sets the entering and leaving variables of :class:`sage.numerical.interactive_simplex_method.LPDictionary` as defined by Dantzig's rule.

        Input:
        - ``current_dictionary`` -- an instance of :class:`sage.numerical.interactive_simplex_method.LPDictionary`

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: pivot = SimplexMethodPivot("dantzig")
        sage: pivot(D)
        sage: D.entering()
        x1
        sage: D.leaving()
        x4
        """
        if current_dictionary.entering() is None:
            entering_index = argmax(current_dictionary.objective_coefficients())
            entering_variable = current_dictionary.nonbasic_variables()[entering_index]
            current_dictionary.enter(entering_variable)
        if current_dictionary.leaving() is None:
            possible = current_dictionary.possible_leaving()
            if possible:
                leaving_variable = min(current_dictionary.ratios())[1]
                current_dictionary.leave(leaving_variable)

class NWRule(AbstractSimplexMethodPivotRule):
    r"""
    Selects the entering and leaving variable based on a Normalzied Weight (NW) pivot rule.

    For a linear program `min(cx: Ax\leq b) = min(cx: x\in P)`, an NW-rule is defined by a normalization `\eta: \mathbb{R}^d\to \mathbb{R}` and a weight `w\in \mathbb{R}^n`.
    Given `v\in P` s.t. `cv \neq cv_{\text{opt}}`, pick the next vertex `u_* = \text{argmax}\{ \frac{w^t(u-v)}{\eta(u-v)}:` `u` adjacent to `v` and `c^tu>c^tv\}`.
    See [BDLS2022]_.

    EXAMPLES::

        sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: pivot = SimplexMethodPivot("NW_rule")
        sage: def eta(v): # use 2 norm as normalization
        ....:     return v.norm()
        sage: pivot(D, False, eta, [1,2])
        sage: D.entering()
        x2
        sage: D.leaving()
        x3

    NW rule is not eqipped to solve both a Phase I and a phase II with :meth:`run_simplex_method`::

        sage: A = ([1, 1], [3, 1], [-1, -1])
        sage: b = (1000, 1500, -400)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: def eta(v):
        ....:     return v.norm()
        sage: P.run_simplex_method('NW_rule', eta, [1,2])
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for *: 'Ambient free module of rank 2 over the principal ideal domain Integer Ring' and 'Vector space of dimension 3 over Rational Field'

    A Phase I and a Prase II problem require different NW rules due to the difference in ambient dimension of the respective problems::

        sage: A = ([1, 1], [3, 1], [-1, -1])
        sage: b = (1000, 1500, -400)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: def eta(v):
        ....:     return v.norm()
        sage: D = P.auxiliary_problem().initial_dictionary() # Phase I
        sage: D.enter('x0')
        sage: D.leave('x5')
        sage: D.update()
        sage: D.run_simplex_method('NW_rule', eta, [1,1,4])
        \begin{equation*}
        ...
        \end{equation*}
        Entering: $x_{2}$. Leaving: $x_{0}$.
        \begin{equation*}
        ...
        \end{equation*}
        sage: d = P.feasible_dictionary(D)
        sage: d.run_simplex_method('NW_rule', eta, [1,1])
        \begin{equation*}
        ...
        \end{equation*}
        Entering: $x_{5}$. Leaving: $x_{3}$.
        \begin{equation*}
        ...
        \end{equation*}
        Entering: $x_{1}$. Leaving: $x_{4}$.
        \begin{equation*}
        ...
        \end{equation*}
    """
    def dictionary_pivot(current_dictionary, eta, weight):
        r"""
        Sets the entering and leaving variables of `:class:`~sage.numerical.interactive_simplex_method.LPDictionary` as defined by Normalized Weight pivot rule.

        Input:
        - ``current_dictionary`` -- an instance of `:class:`~sage.numerical.interactive_simplex_method.LPDictionary`

        - ``eta`` -- callable function from ambient space of ``current_dictionary`` to the real numbers.

        - ``weight`` -- object convertible to a sage ``vector`` in the ambient space of current_dictionary.

        EXAMPLES::

            sage: from sage.numerical.pivot_rules_for_simplex_method import SimplexMethodPivot
            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.entering()

            sage: D.leaving()

            sage: def eta(v): # use 2 norm as normalization
            ....:     return v.norm()
            sage: pivot = SimplexMethodPivot("NW_rule")
            sage: pivot(D, False, eta, [1,1])
            sage: D.entering()
            x1
            sage: D.leaving()
            x4

            Inputs are assumed to be correct for the current problem::
            sage: D.update()
            sage: pivot(D, False, eta, 1)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: eta = [1,2,1]
            sage: pivot(D, False, eta, [1,1])
            Traceback (most recent call last):
            ...
            TypeError: 'list' object is not callable
        """
        eta = eta
        weight = vector(weight)
        if current_dictionary.entering() is None:
            current_vertex = current_dictionary.basic_solution() # This current vertex v of P.
            best_value = InfinityRing(float('-inf'))
            best_vars = None
            for var_in in current_dictionary.possible_entering(): # Explore adjacent vertices to v such that c^tu>c^tv.
                copy_of_dict = deepcopy(current_dictionary) # TODO: custom copy method for LPAbstractDictionary or better way to explore multiple vertices.
                copy_of_dict.enter(var_in)
                for var_out in copy_of_dict.possible_leaving():
                    copy_of_dict.leave(var_out)
                    copy_of_dict.update()
                    next_vertex = copy_of_dict.basic_solution() # Vertex u adjacent to v in gr(P).
                    if eta(next_vertex - current_vertex) != 0:
                        current_value = (weight * (next_vertex - current_vertex)) / eta(next_vertex - current_vertex) # NW-rule
                        if current_value > best_value:
                            best_vars = var_in, var_out
                            best_value = current_value
                    copy_of_dict = deepcopy(current_dictionary)
            if best_vars is not None: # Don't set entering or leaving variables if we haven't found any strict improvement.
                current_dictionary.enter(best_vars[0])
                if current_dictionary.leaving() is None:
                    current_dictionary.leave(best_vars[1])


from sage.numerical.interactive_simplex_method import LPAbstractDictionary
