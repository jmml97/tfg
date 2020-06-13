from sage.coding.bch_code import BCHCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder

DEBUG = false

# General functions

def logger(s):
    if DEBUG:
        print(s)

def _to_complete_list(poly, length):
    r"""
    Returns the vector of length exactly ``length`` corresponding to the
    coefficients of the provided polynomial. If needed, zeros are added.
    INPUT:
    - ``poly`` -- a polynomial
    - ``length`` -- an integer
    OUTPUT:
    - the list of coefficients
    """
    L = poly.coefficients(sparse=False)
    return L + [poly.base_ring().zero()] * (length - len(L))

# Functions for PGZ Decoder

def xentries_m_mu(i, j, s):
    return s[i+j]

def xentries_m_syn(i, j, X):
    return X[j]^(i+1)

class BCHPGZDecoder(Decoder):
    r"""
    Constructs a decoder for BCH Codes based on the PGZ algorithm.
    The decoding algorithm works as follows:
    - First, 
    - Then, 
    - Finally, 
    INPUT:
    - ``code`` -- A code associated to this decoder
    """
    
    def __init__(self, code):
        if not isinstance(code, BCHCode):
            raise ValueError("code has to be a BCHCode")
        self._polynomial_ring = code._polynomial_ring
        self._correction_capability = floor((code.designed_distance() - 1)/2)
        super(BCHPGZDecoder, self).__init__(
            code, code.ambient_space(), "BCHPGZDecoder")
    
    def __eq__(self, other):
        r"""
        Tests equality between BCHPGZDecoder objects.
        """
        return (isinstance(other, BCHPGZDecoder) and
                self.code() == other.code())
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Peterson-Gorenstein-Zierler algorithm based decoder for %s" % self.code()
    
    def decode_to_code(self, word):
        r"""
        Corrects the errors in ``word`` and returns a codeword.
        INPUT:
        - ``r`` -- a codeword of ``self``
        OUTPUT:
        - a vector of ``self``'s message space
        """
        R = self._polynomial_ring
        F = self.code().base_ring()
        logger(R)
        logger(R.gen())
        y = R(word.list())
        t = self.correction_capability()
        g = self.code().generator_polynomial()
        a = self.code().primitive_root()
        
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
            sigma = sigma + sol_mu[len(sol_mu)-1-i]*R.gen()^(i+1)
            
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
        logger('magnitudes de error E: ' + str(E))
        
        # Paso final: hallamos el error y se lo restamos al mensaje recibido
        
        e = R(0)
        for i in [0..mu-1]:
            e = e + E[i]*R.gen()^(k[i])
        
        logger('error e: ' + str(e))
        
        c = y - e
        
        return vector(F, _to_complete_list(c, self.code().length()))
    
    def correction_capability(self):
        return self._correction_capability

BCHCode._registered_decoders["BCHPGZDecoder"] = BCHPGZDecoder