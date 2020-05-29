""" Requiere como argumento un polinomio de la forma x^n - 1
    Devuelve una lista de polinomios generadores de códigos cíclicos
    de longitud n junto a sus conjuntos característicos, indicando la
    raíz primitiva escogida por Sage
"""
def ctos_caracteristicos(pol):
    n = pol.degree()
    gen_cto = []
    for g in generadores(pol):
        C = codes.CyclicCode(generator_pol=g, length=n)
        gen_cto.append((g, C.defining_set(), C.primitive_root()))
    return gen_cto