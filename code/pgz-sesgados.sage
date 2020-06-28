from sage.coding.linear_code import AbstractLinearCode
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

def left_lcm(pols):
    r"""
    Returns the left least common multiple of the polynomials in ``pols``
    INPUT:
    - ``pols`` -- a list of polynomials
    OUTPUT:
    - the left least common multiple of the polynomials
    """
    F = pols[0].base_ring()
    R = pols[0].parent()
    res = R(1)
    for p in pols:
        res = res.left_lcm(p)
    return res

def norm(i, gamma, sigma):
    r"""
    Returns the ``i``th-norm of ``gamma`` with automorphism ``sigma``
    INPUT:
    - ``i`` -- the order of the norm
    - ``gamma`` -- the element to which the norm will be calculated
    - ``sigma`` -- the automorphism used to calculate the norm
    OUTPUT:
    - the ``i``th-norm of ``gamma`` with automorphism ``sigma``
    """
    r = gamma
    if (i == 0):
        return 1
    else:
        for j in (1..i-1):
            r = r*(sigma^j)(gamma)
    return r

# Functions for PGZ algorithm

def xentries_N(i, j, sigma, b):
    return norm(i, (sigma^j)(b), sigma)

def xentries_s_t(i, j, s, sigma, a):
    return (sigma^-j)(s[i+j])*(sigma^i)(a)

def xentries_M_rho(i, j, rho, sigma):
    if (i > j or j > len(rho) - 1 + i):
        return 0
    else:
        return (sigma^i)(rho[j-i])

def xentries_Sigma(i, j, sigma, a, k):
    return (sigma^(k[j]+i))(a)

# Functions for Skew Codes

def xentries_M_g(i, j, g, sigma):
    l = _to_complete_list(g, sigma.order())
    return (sigma^i)(l[j-i])

class SkewCyclicCode(AbstractLinearCode):
    r"""
    Representation of a SkewCyclicCode as a linear code.

    INPUT:
    - ``generator_pol`` -- A generator polynomial of the code
    """
    
    _registered_encoders = {}
    _registered_decoders = {}
    
    def __init__(self, generator_pol=None):
        F = generator_pol.base_ring()
        R = generator_pol.parent()
        n = R.twist_map().order()
        D = generator_pol.degree() + 1
        
        self._polynomial_ring = R
        self._generator_polynomial = generator_pol
        self._primitive_root = F.primitive_element()
        self._length = n
        self._dimension = n - generator_pol.degree()
        
        super(SkewCyclicCode, self).__init__(F, n, "SkewCyclicVectorEncoder", "SkewCyclicDecoder")
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return ("[%s, %s] Skew Cyclic Code on %s"
                % (self.length(), self.dimension(),
                   self.polynomial_ring()))
    
    def generator_polynomial(self):
        r"""
        Returns the generator polynomial of ``self``.
        """
        return self._generator_polynomial
    
    def polynomial_ring(self):
        r"""
        Returns the underlying polynomial ring of ``self``.
        """
        return self._polynomial_ring
    
    def primitive_root(self):
        r"""
        Returns a primitive root of the underlying field of ``self``.
        """
        return self._primitive_root
    
    def ring_automorphism(self):
        r"""
        Returns the ring automorphism of the underlying skew field of ``self``.
        """
        return self._polynomial_ring.twist_map()

class SkewRSCode(SkewCyclicCode):
    r"""
    Representation of a skew RS code.

    We propose two different ways to create a new CyclicCode, either by
    providing:
    - the generator polynomial in a skew polynomial ring
    - a list of b-roots of the underlying field, whose left least common multiple will be the generator polynomial used

    INPUT:
    - ``generator_pol`` -- (default: ``None``) the generator polynomial
      of ``self``. That is, the highest-degree monic polynomial which divides
      every polynomial representation of a codeword in ``self``.
    - ``b_roots`` -- (default: ``None``) a list of b-roots of the underlying field, whose left least common multiple will be the generator polynomial used.
    """
    
    def __init__(self, generator_pol=None, b_roots=None):
        r"""
        TESTS:

        Either if one privides a generator polynomial or we compute it using
        the b-roots, we check that that the generator polynomial divides 
        `x^{n} - 1`, where `n` is the order of the skew polynomial ring's
        automorphism.
        """
        
        if (b_roots is not None and generator_pol is None):
            F = b_roots[0].base_ring()
            R = b_roots[0].parent()
            n = R.twist_map().order()
            
            for b_root in b_roots:
                if not b_root.left_divides(R.gen()^n - 1):
                    raise ValueError("All the b_roots must divide x^n - 1, "
                                 "where n is the automorphism order.")
            
            self._generator_polynomial = left_lcm(b_roots)
            
        elif (b_roots is None and generator_pol is not None):
            F = generator_pol.base_ring()
            R = generator_pol.parent()
            n = R.twist_map().order()
            
            if not generator_pol.left_divides(R.gen()^n - 1):
                    raise ValueError("The generator_pol must divide x^n - 1, "
                                 "where n is the automorphism order.")

            self._generator_polynomial = generator_pol

        else:
            raise AttributeError("You must provide either a list of factors or a"
                                 "generator polynomial.")
        
        self._polynomial_ring = R
        self._primitive_root = F.primitive_element()
        self._length = n
        self._dimension = n - self._generator_polynomial.degree()
        self._designed_distance = self._generator_polynomial.degree() + 1
        
        super(SkewCyclicCode, self).__init__(F, n, "SkewCyclicVectorEncoder", "SkewRSPGZDecoder")
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return ("[%s, %s] Skew RS Code on %s"
                % (self.length(), self.dimension(),
                   self.polynomial_ring()))
    
    def designed_distance(self):
        r"""
        Returns the designed distance of ``self``.
        """
        return self._designed_distance
    
class SkewCyclicVectorEncoder(Encoder):
    r"""
    An encoder which can encode vectors into codewords.
    Let `C` be a skew cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.
    Let `m = (m_1, m_2, \dots, m_k)` be a vector in `F^{k}`.
    To encode `m`, this encoder does the following multiplication:
    `mM(g)`, where M(g) is the generator matrix for g.
    INPUT:
    - ``code`` -- The associated code of this encoder
    """
    
    def __init__(self, code):
        r"""
        TESTS:

        We check that the provided code is a SkewCyclicCode
        """

        if not isinstance(code, SkewCyclicCode):
            raise ValueError("code has to be a SkewCyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super(SkewCyclicVectorEncoder, self).__init__(code)
        
    def __eq__(self, other):
        r"""
        Tests equality between SkewCyclicVectorEncoder objects.
        """
        return (isinstance(other, SkewCyclicVectorEncoder)
            and self.code() == other.code())
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Vector-style encoder for %s" % self.code()
    
    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``
        """
        g = self.code().generator_polynomial()
        sigma = self.code().ring_automorphism()
        n = self.code().length()
        k = n - g.degree()
        return matrix(k, n, lambda i, j: xentries_M_g(i, j, g, sigma))
    
class SkewCyclicPolynomialEncoder(Encoder):
    r"""
    An encoder encoding polynomials into codewords.
    Let `C` be a skew cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.
    This encoder encodes any polynomial `p \in F[x]_{<k}` by computing
    `c = p \times g` and returning the vector of its coefficients.
    INPUT:
    - ``code`` -- The associated code of this encoder
    """
    
    def __init__(self, code):
        r"""
        TESTS:

        We check that the provided code is a SkewCyclicCode
        """

        if not isinstance(code, SkewCyclicCode):
            raise ValueError("code has to be a SkewCyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super(SkewCyclicPolynomialEncoder, self).__init__(code)
        
    def __eq__(self, other):
        r"""
        Tests equality between SkewCyclicPolynomialEncoder objects.
        """
        return (isinstance(other, SkewCyclicPolynomialEncoder)
            and self.code() == other.code())
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Polynomial-style encoder for %s" % self.code()
    
    def message_space(self):
        r"""
        Returns the message space of ``self``.
        """
        return self._polynomial_ring
    
    def encode(self, p):
        r"""
        Transforms ``p`` into an element of the associated code of ``self``.
        INPUT:
        - ``p`` -- A polynomial from ``self`` message space
        OUTPUT:
        - A codeword in associated code of ``self``
        """
        C = self.code()
        k = C.dimension()
        n = C.length()
        if p.degree() >= k:
            raise ValueError("Degree of the message must be at most %s" % k - 1)
        res = _to_complete_list(p * C.generator_polynomial(), n)
        return vector(C.base_field(), res)
    
    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.
        INPUT:
        - ``c`` -- A vector with the same length as the code
        OUTPUT:
        - An element of the message space
        """
        R = self.message_space()
        g = self.code().generator_polynomial()
        p = R(c.list())
        return p // g

class SkewCyclicDecoder(Decoder):
    r"""
    Constructs a decoder for Skew Cyclic Codes.
    The decoding algorithm IS NOT IMPLEMENTED.
    INPUT:
    - ``code`` -- A code associated to this decoder
    """
    
    def __init__(self, code):
        r"""
        TESTS:

        We check that the provided code is a SkewCyclicCode
        """

        if not isinstance(code, SkewCyclicCode):
            raise ValueError("code has to be a SkewCyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super(SkewCyclicDecoder, self).__init__(
            code, code.ambient_space(), "SkewCyclicDecoder")
    
    def __eq__(self, other):
        r"""
        Tests equality between SkewCyclicDecoder objects.
        """
        return (isinstance(other, SkewCyclicDecoder) and
                self.code() == other.code())
    
    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Decoder for %s" % self.code()
    
    def decode_to_code(self, word):
        r"""
        Corrects the errors in ``word`` and returns a codeword.
        INPUT:
        - ``r`` -- a codeword of ``self``
        OUTPUT:
        - a vector of ``self``'s message space

        THIS METHOD IS NOT IMPLEMENTED
        """
        return vector(0)
    
class SkewRSPGZDecoder(Decoder):
    r"""
    Constructs a decoder for Skew RS Codes based on the PGZ algorithm for Skew Cyclic codes.
    The decoding algorithm works as follows:
    - First, the syndromes are computed
    - Then, the error locator polynomial is computed using the syndromes, 
      obtaining the error locations
    - Finally, the syndrome's equation system is solved, obtaining the error 
      and the decoded codeword
    INPUT:
    - ``code`` -- A code associated to this decoder
    """
    
    def __init__(self, code):
        if not isinstance(code, SkewRSCode):
            raise ValueError("code has to be a SkewRSCode")
        self._polynomial_ring = code._polynomial_ring
        self._correction_capability = floor((code.designed_distance() - 1)/2)
        super(SkewRSPGZDecoder, self).__init__(
            code, code.ambient_space(), "SkewRSPGZDecoder")
    
    def __eq__(self, other):
        r"""
        Tests equality between SkewRSPGZDecoder objects.
        """
        return (isinstance(other, SkewRSPGZDecoder) and
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
        y = R(word.list())
        y_list = _to_complete_list(y, self.code().length())
        sigma = self.code().ring_automorphism()
        a = self.code().primitive_root()
        b = sigma(a)*a^-1
        t = self.correction_capability()
        n = self.code().length()

        N = matrix(n, lambda i, j: xentries_N(i, j, sigma, b))

        # Step 1: compute all the syndromes
        s = []
        i = 0
        for i in (0..2*t - 1):
            s_i = 0
            for j in (0..n-1):
                s_i = s_i + y_list[j]*norm(j, (sigma^i)(b), sigma)
            s.append(s_i)
        logger("s, syndromes vector: ", s)

        # Check if syndromes are all zero
        if all(s_i == 0 for s_i in s):
            return 0

        # Step 2: compute the error locator polynomial and the error locations
        # Compute the matrix S_t
        S_t = matrix(t+1, t, lambda i, j: xentries_s_t(i, j, s, sigma, a))
        logger("S_t:\n", S_t)

        rcef_S_t = S_t.transpose().rref().transpose()
        logger("rcef_S_t:\n", rcef_S_t)

        mu = 0
        for i in (1..rcef_S_t.ncols()):
            sub = rcef_S_t[list(range(i)),list(range(i))]
            if sub == matrix.identity(F, i):
                mu = i

        # Obtain the possible error locator polynomial from the a_i
        # coefficients from the matrix rcef_S_t
        rho = [-a_i for a_i in rcef_S_t[[mu], list(range(0, mu))].list()]
        rho.append(1)
        logger("rho: ", rho)

        rho_ = _to_complete_list(R(rho), n)
        rho_N = vector(rho_)*N
        logger("rho_N: ", rho_N)

        # Compute the error locations
        k = [i for i, e in enumerate(rho_N) if e == 0]
        v = len(k)

        # If we didn't find all the error locations
        if mu != v:
            logger("Case mu != v")
            M_rho = matrix(n - mu, n, lambda i, j: xentries_M_rho(i, j, rho, sigma))
            logger("M_rho:\n", M_rho)
            N_rho = M_rho*N
            logger("N_rho:\n", N_rho)
            H_rho = N_rho.rref()
            logger("H_rho:\n", H_rho)

            H_ = H_rho

            # Rows from H_rho that are not a canonical vector epsilon_i are 
            # deleted
            for i in range(H_rho.nrows()):
                sub = H_rho.row(i)
                epsilon = matrix.identity(F, H_rho.ncols()).rows()
                if (sub not in epsilon):
                    H_ = H_.delete_rows([i])
            logger("H_:\n", H_)

            # New error locations are computed
            k = [i for i, e in enumerate(H_.columns()) if e == vector([0] * H_.nrows())]
            v = len(k)

        logger("k: ", k)
        logger("v: ", v)
        logger("Note: solve for E, where E*Sigma.transpose() = b_syn")
        Sigma = matrix(v, lambda i, j: xentries_Sigma(i, j, sigma, a, k))
        logger("Sigma:\n", Sigma)

        b_syn = [a*s[0]]
        for i in (1..v-1):
            b_syn.append(sigma(a)*s[i])

        logger("b_syn: ", b_syn)

        # Step 3: solve the syndrome's system
        E = Sigma.transpose().solve_left(vector(b_syn))
        logger("E: ", E)

        # Step 4: we compute the error and substract it from the received 
        # message, obtaining the decoded codeword
        e = R(0)
        for i in range(len(k)):
            e = e + E[i]*R.gen()^(k[i])
        logger("error: ", e)
        logger("m = y - e: ", y - e)

        return vector(_to_complete_list(y - e, self.code().length()))
    
    def correction_capability(self):
        r"""
        Returns the correction capability of ``self``.
        """
        return self._correction_capability
    
SkewCyclicCode._registered_encoders["SkewCyclicVectorEncoder"] = SkewCyclicVectorEncoder
SkewCyclicCode._registered_encoders["SkewCyclicPolynomialEncoder"] = SkewCyclicPolynomialEncoder
SkewCyclicCode._registered_decoders["SkewCyclicDecoder"] = SkewCyclicDecoder

SkewRSCode._registered_decoders["SkewRSPGZDecoder"] = SkewRSPGZDecoder