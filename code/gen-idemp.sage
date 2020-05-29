from itertools import chain, combinations

# -- https://stackoverflow.com/a/52225720/1629814
""" Nos da todas las posibles combinaciones de elementos de una lista
    Ejemplo: 
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
"""
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

"""Multiplica todos los elementos del iterable y devuelve el resultado"""
def mult(iterable):
    rv = 1
    for k in iterable:
        rv = rv * k
    return rv
# -- 

""" Requiere como argumento un polinomio de la forma x^n - 1
    Devuelve una lista de polinomios generadores de códigos cíclicos
    de longitud n
"""
def generadores(pol):
    factores = list(pol.factor())
    factores = [factor[0] for factor in factores]
    k = powerset(factores)
    result = [1]
    result.extend(map(mult, (x for x in k if len(x)>0)))
    return result 

""" Requiere como argumento un polinomio de la forma x^n - 1
    Devuelve una lista de tuplas (g, e), donde g es un polinomio generador
    de un código cíclico de longitud n y e es el idempotente correspondiente
"""

def generadores_idempotentes(pol):
    gen = generadores(pol)
    gen_idempotentes = []
    
    for g in gen:
        d, s, t = xgcd(g, pol/g)
        gen_idempotentes.append((g, s*g))
    
    return gen_idempotentes