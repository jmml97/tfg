# *****************************************************************************
#             Copyright (C) 2020 José María Martín Luque
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************
from itertools import chain, combinations

# -- https://stackoverflow.com/a/52225720/1629814

def powerset(iterable):
    r"""
    Returns all possible combinations of elements in ``iterable``.
    INPUT:
    - ``iterable`` -- An iterable container
    OUTPUT:
    - All possible combinations of elements in ``iterable``
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def mult(iterable):
    r"""
    Returns the product of all elements in ``iterable``.
    INPUT:
    - ``iterable`` -- An iterable container
    OUTPUT:
    - The product of all elements in ``iterable``
    """
    rv = 1
    for k in iterable:
        rv = rv * k
    return rv
# -- 

def generators(poly):
    r"""
    Returns a list of generator polynomials of cyclic codes of length n, where 
    n is the degree of ``poly``.
    INPUT:
    - ``poly`` -- A polynomial of the form x^n - 1
    OUTPUT:
    - A list of generator polynomials of cyclic codes of length n
    """
    factors = list(poly.factor())
    factors = [factor[0] for factor in factors]
    k = powerset(factors)
    result = [1]
    result.extend(map(mult, (x for x in k if len(x)>0)))
    return result 

def generator_and_idempotents(poly):
    r"""
    Returns a list of generator polynomials of cyclic codes of length n, where 
    n is the degree of ``poly``, with their corresponding generator idempotent.
    INPUT:
    - ``poly`` -- A polynomial of the form x^n - 1
    OUTPUT:
    - A list of tuples containing a generator polynomial, and a generator idempotent
    """
    gen = generators(poly)
    gen_and_idemp = []
    
    for g in gen:
        d, s, t = xgcd(g, poly/g)
        gen_and_idemp.append((g, s*g))
    
    return gen_and_idemp

def defining_sets(poly):
    r"""
    Returns a list of generator polynomials of cyclic codes of length n, where 
    n is the degree of ``poly``, with their corresponding defining sets and primitive roots.
    INPUT:
    - ``poly`` -- A polynomial of the form x^n - 1
    OUTPUT:
    - A list of tuples containing a generator polynomial, a defining set and a 
    primitive root
    """
    n = poly.degree()
    generators_and_sets = []
    for g in generators(poly):
        C = codes.CyclicCode(generator_pol=g, length=n)
        generators_and_sets.append((g, C.defining_set(), C.primitive_root()))
    return generators_and_sets