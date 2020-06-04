DEBUG = false

def logger(s):
    if DEBUG:
        print(s)

def xentries_m_mu(i, j, s):
    return s[i+j]

def xentries_m_syn(i, j, X):
    return X[j]^(i+1)

def PGZ(C, y):
    delta = C.designed_distance()
    t = floor((delta - 1)/2)
    g = C.generator_polynomial()
    a = C.primitive_root()
    
    logger('polinomio generador: ' + str(g))
    logger('raíz primitiva: ' + str(a))
    
    # Paso 1: calcular los síndromes
    s = []
    for i in range(0, 2*t):
        s.append(y(a^(i+1)))
    
    logger('síndromes: ' + str(s))
        
    # Paso 2: probar mu partiendo de t hasta encontrar matriz m_mu
    # que sea no singular. Resolver el sistema para encontrar sigma(x)
    mu = t
    
    m_mu = matrix(mu, lambda i, j: xentries_m_mu(i, j, s))
    
    while m_mu.determinant() == 0:
        mu = mu - 1
        m_mu = matrix(mu, lambda i, j: xentries_m_mu(i, j, s))
    
    logger('tamaño de m_mu: ' + str(mu))
    logger('matriz m_mu: \n' + str(m_mu))
    
    if m_mu == matrix():
        print("Se produjeron más de ", t, " errores: no podemos decodificar")
        return
    
    b_mu = []
    
    for i in [mu..2*mu-1]:
        b_mu.append(-s[i])
        
    b_mu = vector(b_mu)
    
    logger('vector b_mu: ' + str(b_mu))
    
    sol_mu = m_mu.solve_right(b_mu)
    
    logger('matriz de soluciones de m_mu*S = b_mu: ' + str(sol_mu))
    
    sigma = 1
    l = len(sol_mu) - 1
    for i in [0..l]:
        sigma = sigma + sol_mu[len(sol_mu)-1-i]*x^(i+1)
        
    logger('polinomio localizador sigma(x): ' + str(sigma))
    
    # Paso 3: encontrar las raíces de sigma(x) e invertirlas para
    # encontrar los números de posición de error X_j
    # Obtenemos también los k_j
    r = sigma.roots()
    
    logger('raíces de sigma(x): ' + str(r))
    
    if r == []:
        print('sigma(x) no tiene raíces, no podemos seguir')
        return
    
    l = len(r) - 1
    X = []
    for i in [0..l]:
        X.append(r[i][0]^-1)
        
    logger('X_j: ' + str(X))
        
    k = []
    for i in [0..l]:
        k.append(log(X[i], a))
    
    logger('k_j: ' + str(k))
        
    # Paso 4: Resolvemos las primeras mu ecuaciones del sistema de los
    # síndromes para hallar las magnitudes de error E_j
    
    m_syn = matrix(mu, lambda i, j: xentries_m_syn(i, j, X))
    
    b_syn = []
    for i in [0..mu-1]:
        b_syn.append(s[i])
    b_syn = vector(b_syn)
    
    E = m_syn.solve_right(b_syn)
    
    # Paso final: hallamos el error y se lo restamos al mensaje recibido
    
    e = 0
    for i in [0..mu-1]:
        e = e + E[i]*x^(k[i])
    
    logger('error e: ' + str(e))
    
    c = y - e
    
    return c