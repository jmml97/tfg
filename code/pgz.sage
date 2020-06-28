from sage.coding.bch_code import BCHCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder

DEBUG = false

# General functions

def logger(*s):
    if DEBUG:
        print(*s, sep = '')

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
    - First, syndromes are computed
    - Then, the error locator polynomial is computed
    - Finally, the roots of the error locator polynomial are used to determine the error, which then is used to decode the message
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
        y = R(word.list())
        t = self.correction_capability()
        g = self.code().generator_polynomial()
        a = self.code().primitive_root()
        
        logger('generator polynomial: ', str(g))
        logger('primitive root: ', str(a))
        
        # Step 1: compute syndromes
        s = []
        for i in range(0, 2*t):
            s.append(y(a^(i+1)))
        
        logger('syndromes: ', str(s))
        
        # Step 2: begin trial and error process to find m_mu, beginning with
        # mu = t, until m_mu is non-singular.
        # Then solve, the system to find sigma(x)
        mu = t
        
        m_mu = matrix(mu, lambda i, j: xentries_m_mu(i, j, s))
        
        while m_mu.determinant() == 0:
            mu = mu - 1
            m_mu = matrix(mu, lambda i, j: xentries_m_mu(i, j, s))
        
        logger('m_mu size: ', str(mu))
        logger('matrix m_mu: \n', str(m_mu))
        
        if m_mu == matrix():
            print("More than ", t, "errors have occured in transmission: it's not possible to decode")
            return
        
        b_mu = []
        
        for i in [mu..2*mu-1]:
            b_mu.append(-s[i])
            
        b_mu = vector(b_mu)
        
        logger('vector b_mu: ', str(b_mu))
        
        sol_mu = m_mu.solve_right(b_mu)
        
        logger('solutions of m_mu*S = b_mu: ', str(sol_mu))
        
        sigma = 1
        l = len(sol_mu) - 1
        for i in [0..l]:
            sigma = sigma + sol_mu[len(sol_mu)-1-i]*R.gen()^(i+1)
            
        logger('error locator polynomial sigma(x): ', str(sigma))
        
        # Step 3: find the roots of sigma(x) and invert them to find the
        # error location numbers X_j
        # Calculate the error locations k_j
        r = sigma.roots()
        
        logger('sigma(x) roots: ' + str(r))
        
        if r == []:
            return('sigma(x) has no roots: unable to continue')
        
        l = len(r) - 1
        X = []
        for i in [0..l]:
            X.append(r[i][0]^-1)
            
        logger('X_j: ', str(X))
            
        k = []
        for i in [0..l]:
            k.append(log(X[i], a))
        
        logger('k_j: ', str(k))
        
        # Step 4: solve the first mu equations of the syndromes' system to find error magnitudes E_j
        
        m_syn = matrix(mu, lambda i, j: xentries_m_syn(i, j, X))
        
        b_syn = []
        for i in [0..mu-1]:
            b_syn.append(s[i])
        b_syn = vector(b_syn)
        
        E = m_syn.solve_right(b_syn)
        logger('error magnitudes E: ', str(E))
        
        # Final step: error is calculated and substracted from the received message to find the codeword
        
        e = R(0)
        for i in [0..mu-1]:
            e = e + E[i]*R.gen()^(k[i])
        
        logger('error e: ' + str(e))
        
        c = y - e
        
        return vector(F, _to_complete_list(c, self.code().length()))
    
    def correction_capability(self):
        r"""
        Returns the correction capability of the algorithm for the code associated to `self`.
        """
        return self._correction_capability

BCHCode._registered_decoders["BCHPGZDecoder"] = BCHPGZDecoder