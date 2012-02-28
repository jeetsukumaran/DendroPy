#! /usr/bin/env python

"""
Matrix and Vector Operations.

Original Source:

    Basic Table, Matrix and Vector functions for Python 2.2
    License: Public Domain
    Author: Raymond Hettinger
    email: python@rcn.com
    Updates and documentation:  http://users.rcn.com/python/download/python.htm

Adapted, extended and modified for inclusion in the DendroPy library.

A 100% pure Python module for vector, matrix, and table math operations which
runs on Python 2.1 or later. It neither conflicts with nor requires NumPy.

For matrices, it includes least matrix multiplication, division, least squares
solutions, the LU and QR factorizations, inverses, determinants, and eigenvalue
determination (for real matrices). Vector operations include dot product, outer
product, cross product, L2 norm, and polynomial evaluation. Table operations
provide implicit looping and convenient expression of functions applied to
whole tables. For ease of use, these operations provide broadcasting for
convenient expression of vector/scalar, matrix/scalar, and matrix/vector
operations. In addition, several helper functions are provided to create
elementary matrices and for polynomial and rational curve-fitting.

Organized into a class Table, two sub-classes Matrix and Vec, and assorted
helper functions.

REFERENCES

    [Golub] Gene H. Golub and Charles F. Van Loan (1996) Matrix Computations 3rd Edition.
    [Kincaid] David Kincaid and Ward Cheney (1996) Numerical Recipes 2nd Edition.
    [Mathews] John H. Mathews and Kurtis D. Fink (1999) Numerical Methods Using Matlab 3rd Edition.
    [Press] William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery (1993) Numerical Recipes in C: The Art of Scientific Computing 2nd Edition.

ALGORITHMS

    - Determinants are computed as the product of diagonal elements on the
      upper-triangular matrix in a QR decomposition. A=QR, so det A=det Q * det
      R, and det Q=1, so det A=det R. When R is upper-triangular, det R=prod(
      diag( R ) ).

    - Solve applies a custom solver to a matrix and performs iterative
      improvement [Press section 2.5] up to 10 times or until the solution
      stops improving. Solving Ax=b, and improve to x + solve(A, b-Ax). The
      custom solver is back-substitution for an upper triangular matrix ,
      forward-substitution for a lower triangular matrix, and QR decomposition
      for general matrices. Ax=b is transformed to QRx=b then Rx=Qt b which is
      then solved by back-substitution. This method works for non-square
      matrices and returns the least squares solution when Ax=b is inconsistent
      or over-determined. Note, [Golub 3.5.3] comments that this implementation
      of iterative improvement is unlikely to increase the accuracy of the
      solution since there are few if any significant digits in b-Ax; however,
      my experiments do show improvement. To guard against loss of precision,
      the last step that does not produce an improvement (as measured by the
      norm squared) is thrown away.

    - The QR decomposition uses Householder reflections [Golub algorithm 5.1.1]
      applied in vector form rather than matrix form [Golub section 5.1.4].

    - The eigs routine first makes a similarity transformation (keeping
      eigenvalues constant) to upper Hessenberg form [Golub algorithm
       7.4.2] and then applies the QR method with shifts to force
      an eigenvalue into the lower right corner. The matrix is then
      deflated by one row and column and the method is repeated until all
      eigenvalues are found [Kincaid section 5.5].

    - Function approximation defaults to Chebyshev nodes rather than equally
      spaced points in order to avoid the Runge phenomenon [Matthews section
      4.5] where approximation accuracy decreases when the number of
      sample points are increased. The nodes are located at t(n)=cos(
      (2N+1-2k)*Pi/(2N+2) ) for k=0, 1, ... N. The nodes are then
      rescaled from [-1,1] to the specified interval.

    - Polynomial evaluation uses Horner's method of nested multiplication
      [Kincaid section 1.2].

LIMITATIONS

    - The determinant has a 50% chance of having the wrong sign because the QR
      decomposition may be a reflection of R rather than just a rotation.

    - Matrix sizes are pre-computed when the matrix is constructed. That size
      does not update if an element is deleted. Recommend building a new matrix
      from the desired vectors.

    - The LU decomposition does not use pivoting or scaling resulting in less
      numerical stability and inability to bypass zero elements on the
      diagonal.

    - The QR decomposition results in a Q that is slightly off from being
      orthonormal because it prioritizes the other post-conditions (that A=Q*R
      and R is upper-triangular). As a result the post-condition
      assertion sometimes fails for QtQ=I within a tight tolerance. Either
      lower the tolerance or disable post-condition assertion checking.

    - The QR decomposition has not been tested for matrices with complex
      values.

    - The .eigs() method handles only real matrices with real eigenvalues.

    - Matrix exponentation is only defined for positive integer powers.

"""

import operator, math, random
NPRE, NPOST = 0, 0                    # Disables pre and post condition checks

def pooled_covariance(u, v, population_variance=False):
    """
    Returns pooled covariance matrix of u, v.
    """
    assert len(u[0]) == len(v[0]), "Number of columns in matrices not equal"
    nrow1 = len(u)
    nrow2 = len(v)
    total_rows = nrow1 + nrow2
    f1 = float(nrow1) / total_rows
    f2 = float(nrow2) / total_rows
    s1 = u.covariance_by_cols(population_variance=population_variance)
    s2 = v.covariance_by_cols(population_variance=population_variance)
    pooled_cov = []
    for i, r1 in enumerate(s1):
        pooled_cov.append([])
        for j, c1 in enumerate(r1):
            pooled_cov[-1].append(f1 * s1[i][j] + f2 *s2[i][j])
    pooled_cov = new_matrix(pooled_cov)
    return pooled_cov


def iszero(z):
    """This predicate evaluates to true when x is nearly zero. It is used
    equality testing in the post-condition assertion checks and for determining
    when a row has been zeroed in the eigenvalue routine."""
    return abs(z) < .000001

def getreal(z):
    try:
        return z.real
    except AttributeError:
        return z

def getimag(z):
    try:
        return z.imag
    except AttributeError:
        return 0

def getconj(z):
    try:
        return z.conjugate()
    except AttributeError:
        return z

separator = [ '', '\t', '\n', '\n----------\n', '\n===========\n' ]

class Table(list):
    """
    Table extends UserList to include math operations which apply to every
    element with implicit looping while retaining other useful list features
    such as slices. It is best described by example::

        >>> inventoryStart=Table( [100, 80, 90, 200, 105] ) # Creates object from a list
        >>> purchases=Table( [75, 35, 23, 41, 30] )
        >>> sales=Table( [19,56, 22,60,65] )
        >>> inventoryEnd=inventoryStart + purchases - sales # Note implicit looping
        >>> inventoryEnd
        [156, 59, 91, 181, 70]
        >>> markup=1.25
        >>> unitPrice=unitCost * markup # Note scalar broadcast to list
        >>> unitPrice
        [1.375, 1.3125, 1.4749999999999999, 2.0, 1.1500000000000001]
        >>> totalPurchases=(unitCost * purchases).sum() # Sum of prices times quantities
        >>> totalPurchases
        239.59
        >>> totalSales=(unitPrice * sales).sum()
        >>> totalSales
        326.8249999999999

    To create a Table object, the Table constructor with a list argument. For
    two-dimensional arrays, call the constructor with of list of one
    dimensional Table objects. Example::

        >>> a=Table( [1,2,3,4] )
        >>> b=Table( [2,3,5,7] )
        >>> c=a * 3
        >>> arr=Table( [a,b,c] ) # 2-D array consisting of three 1-D arrays
        >>> arr
        [[1, 2, 3, 4], [2, 3, 5, 7], [3, 6, 9, 12]]
        >>> print a.dim, arr.dim # show array dimensions
        1 2
        >>> arr[2][3] # access elements using [row][col] references
        12
        >>> arr.prod() # Compute the product of all elements
        9797760.0
        >>> arr - Table( [c,b,a] ) # Subtract one 2-D array from another
        [[-2, -4, -6, -8], [0, 0, 0, 0], [2, 4, 6, 8]]
        >>> arr * 3 # Multiply all elements by a scalar
        [[3, 6, 9, 12], [6, 9, 15, 21], [9, 18, 27, 36]]
        >>> arr + b # Broadcast a row to every row of a 2-D array
        [[3, 5, 8, 11], [4, 6, 10, 14], [5, 9, 14, 19]]


    As the addition operator ('+') has been overridden to perform table addition,
    the method `concat` is provided to make a new table object with the new element
    added to the end::

        >>> [1,2,3,4] + [5] # With lists, the plus operator appends to a list
        [1, 2, 3, 4, 5]
        >>> Table([1,2,3,4]) + 5 # Table object use plus for addition
        [6, 7, 8, 9]
        >>> Table([1,2,3,4]).concat( [5] ) # So a concat method is needed for appending
        [1, 2, 3, 4, 5]

    """
    dim = 1
    concat = list.__add__      # A substitute for the overridden __add__ method

    def __getslice__( self, i, j ):
        return self.__class__( list.__getslice__(self,i,j) )

    def __init__( self, elems ):
        list.__init__( self, elems )
        if len(elems) and hasattr(elems[0], 'dim'): self.dim = elems[0].dim + 1

    def __str__( self ):
        return separator[self.dim].join( map(str, self) )

    def map( self, op, rhs=None ):
        """
        Apply a unary operator to every element in the matrix or a binary operator to corresponding
        elements in two arrays.  If the dimensions are different, broadcast the smaller dimension over
        the larger (i.e. match a scalar to every element in a vector or a vector to a matrix).

        (Unary operator) Apply a function of one variable to every element::

            >>> b.map( math.cos )
            [-0.4161468365471, -0.989992496600, 0.283662185463, 0.75390225434]
            >>> arr.map( lambda x: x%2 and 3*x+1 or x/2 )
            [[4, 1, 10, 2], [1, 10, 16, 22], [10, 3, 28, 6]]

        (Binary operator) Apply a function of two variables to corresponding
        elements in two arrays. If the arrays are not of the same dimension,
        broadcasting occurs (rows or scalars are repeated to match the array of
        higher dimension).

            >>> a.map( lambda x,y: 10*x+y, b ) # Pair elements of two 1-D arrays
            [12, 23, 35, 47]
            >>> a.map( lambda x,y: 10*x+y, 5 ) # Broadcast the 5 to every element of a 1-D array
            [15, 25, 35, 45]
            >>> arr.map( lambda x,y: 10*x+y, b ) # Broadcast a 1D row to every row of a 2-D array
            [[12, 23, 35, 47], [22, 33, 55, 77], [32, 63, 95, 127]]

        """
        if rhs is None:                                                 # Unary case
            return self.dim==1 and self.__class__( map(op, self) ) or self.__class__( [elem.map(op) for elem in self] )
        elif not hasattr(rhs,'dim'):                                    # List / Scalar op
            return self.__class__( [op(e,rhs) for e in self] )
        elif self.dim == rhs.dim:                                       # Same level Vector / Vector or Matrix / Matrix
            assert NPRE or len(self) == len(rhs), 'Table operation requires len sizes to agree'
            return self.__class__( map(op, self, rhs) )
        elif self.dim < rhs.dim:                                        # Vector / Matrix
            return self.__class__( [op(self,e) for e in rhs]  )
        return self.__class__( [op(e,rhs) for e in self] )         # Matrix / Vector

    def __mul__( self, rhs ):
        return self.map( operator.mul, rhs )

    def __div__( self, rhs ):
        return self.map( operator.div, rhs )

    def __sub__( self, rhs ):
        return self.map( operator.sub, rhs )

    def __add__( self, rhs ):
        return self.map( operator.add, rhs )

    def __rmul__( self, lhs ):
        return self*lhs

    def __rdiv__( self, lhs ):
        return self*(1.0/lhs)

    def __rsub__( self, lhs ):
        return -(self-lhs)

    def __radd__( self, lhs ):
        return self+lhs

    def __abs__( self ):
        return self.map( abs )

    def __neg__( self ):
        return self.map( operator.neg )

    def conjugate( self ):
        return self.map( getconj )

    def real( self ):
        return self.map( getreal  )

    def imag( self ):
        return self.map( getimag )

    def flatten( self ):
        """
        Collapse an table object to a 1D list.

            >>> arr.flatten()
            [3, 6, 9, 12, 2, 3, 5, 7, 1, 2, 3, 4]

        """
        if self.dim == 1:
            return self
        return reduce( lambda cum, e: e.flatten().concat(cum), self, [] )

    def prod( self ):
        return reduce(operator.mul, self.flatten(), 1.0)

    def sum( self ):
        return reduce(operator.add, self.flatten(), 0.0)

    def exists( self, predicate ):
        for elem in self.flatten():
            if predicate(elem):
                return 1
        return 0

    def forall( self, predicate ):
        """
        Checks every element to see if all instances are true. The predicate is a boolean function of one variable.

            >>> def isprime(x): return x==2 or 2**(x-1)%x==1
            ...
            >>> a.forall( isprime )
            0
            >>> b.forall(isprime)
            1

        """
        for elem in self.flatten():
            if not predicate(elem):
                return 0
        return 1

    def __eq__( self, rhs ):
        return (self - rhs).forall( iszero )

class Vector(Table):
    """
    The Vector class extends Table and is intended for 1-D arrays. It adds pretty
    printing and standard math vector operations. It is best described by
    example.

        >>> from matfunc import Vec
        >>> a=Vec( [4,-2, 5] ) # Constructed just like the Table class
        >>> b=Vec( [3,10,-6] )
        >>> a
        [4, -2, 5] # Normal table representation
        >>> print a # Pretty printed with built-in __str__() method
        4.000 -2.000 5.000
        >>> a.dot(b) # Computes the dot or inner product of 'a' and 'b'
        -38.0
        >>> a.norm() # Computs the length, hypoteneuse or L2 norm of 'a'
        6.7082039324993694
        >>> print a.normalize() # Vector of length 1 in same direction as 'a'
        0.596 -0.298 0.745
        >>> print a.outer(b) # Constructs a multiplication table from 'a' and 'b'
        12.000 40.000 -24.000
        -6.000 -20.000 12.000
        15.000 50.000 -30.000
        >>> print a.cross(b) # The vector cross product in perpendicular to 'a' and 'b'
        -38.000 39.000 46.000
        >>> Vec([6,3,4]).polyval(5) #evaluates 6*x**2 + 3*x + 4 at x=5
        169.0

    The .house(index) method is a helper for the matrix class and is used to
    compute Householder reflection vectors which are used to zero out parts of
    a vector.
    """

    def dot( self, otherVector ):
        return reduce(operator.add, map(operator.mul, self, otherVector), 0.0)

    def norm( self ):
        return math.sqrt(abs( self.dot(self.conjugate()) ))

    def normalize( self ):
        return self / self.norm()

    def outer( self, otherVector ):
        return new_matrix([otherVector*x for x in self])

    def cross( self, otherVector ):
        'Compute a Vector or Cross Product with another vector'
        assert len(self) == len(otherVector) == 3, 'Cross product only defined for 3-D vectors'
        u, v = self, otherVector
        return Vector([ u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] ])

    def house( self, index ):
        'Compute a Householder vector which zeroes all but the index element after a reflection'
        v = Vector( Table([0]*index).concat(self[index:]) ).normalize()
        t = v[index]
        sigma = 1.0 - t**2
        if sigma != 0.0:
            t = v[index] = t<=0 and t-1.0 or -sigma / (t + 1.0)
            v /= t
        return v, 2.0 * t**2 / (sigma + t**2)

    def polyval( self, x ):
        'Vector([6,3,4]).polyval(5) evaluates to 6*x**2 + 3*x + 4 at x=5'
        return reduce( lambda cum,c: cum*x+c, self, 0.0 )

    def ratval( self, x ):
        'Vector([10,20,30,40,50]).ratfit(5) evaluates to (10*x**2 + 20*x + 30) / (40*x**2 + 50*x + 1) at x=5.'
        degree = len(self) / 2
        num, den = self[:degree+1], self[degree+1:] + [1]
        return num.polyval(x) / den.polyval(x)

class Matrix(Table):
    """
    Matrix routines are organized into a class hierarchy starting with Matrix, its
    subclass `Square` which has a `Triangular` subclass which splits into `UpperTri` and
    `LowerTri`. This organization assigns efficient algorithms to each structure
    (i.e. solving upper triangular matrices with back-substitution) and adds
    structure specific methods (i.e. eigenvalues and inverses only have meaning for
    square matrices).

    `new_matrix(listOflists, type)` is a factory method for creating new
    matrices from a list of lists or a list of `Vector` objects. It
    automatically selects the `Matrix` class for rectangular matrices, `Square`
    when the row and column lengths agree, and `LowerTri` or `UpperTri` if
    there are blocks of zeroes.

    The methods `tr` (for transpose), `star` (for Hermetian adjoints), `diag`,
    `trace`, and augment all have their usual mathematical meanings. Matrix
    multiplication is accomplished by `mmul` and matrix division by `solve(b)`.
    When the systems are inconsistent or overdetermined, the least squares
    solution is computed.

    `qr` returns the QR decompostion where Q is orthonormal (a rotation in space
    where QtQ=I) and R is upper triangular. If the original matrix A is
    m-by-n, then A=Q*R, Q is m-by-n, and R is n-by-n.

    Square matrices add several other methods. `det` and `inverse` have their normal
    mathematical meanings. `hessenberg` returns a similar matrix (having the same
    eigenvalues) where all elements below the first sub-diagonal are zero.
    `eigs` estimates the eigenvalues of real matrices.

    `lu` returns the LU decomposition of square matrices where L is lower triangular,
    U is upper triangular, and L*U is the original matrix.

    `__pow__` raises matrices to positive integer powers; for example, A**3 gives
    A*A*A. Note, this may be unexpected because all other math operators are
    applied table; for example, A+3 adds 3 to every element individually.
    """
    __slots__ = ['size', 'rows', 'cols']

    def __init__( self, elems ):
        'Form a matrix from a list of lists or a list of Vectors'
        Table.__init__( self, hasattr(elems[0], 'dot') and elems or map(Vector,map(tuple,elems)) )
        self.size = self.rows, self.cols = len(elems), len(elems[0])

    def tr( self ):
        'Tranpose elements so that Transposed[i][j] = Original[j][i]'
        return new_matrix(zip(*self))

    def star( self ):
        'Return the Hermetian adjoint so that Star[i][j] = Original[j][i].conjugate()'
        return self.tr().conjugate()

    def diag( self ):
        'Return a vector composed of elements on the matrix diagonal'
        return Vector( [self[i][i] for i in range(min(self.size))] )

    def trace( self ):
        return self.diag().sum()

    def mmul( self, other ):
        'Matrix multiply by another matrix or a column vector '
        if other.dim==2:
            return new_matrix( map(self.mmul, other.tr()) ).tr()
        assert NPRE or self.cols == len(other)
        return Vector( map(other.dot, self) )

    def augment( self, otherMat ):
        'Make a new matrix with the two original matrices laid side by side'
        assert self.rows == otherMat.rows, 'Size mismatch: %s * %s' % (`self.size`, `otherMat.size`)
        return new_matrix( map(Table.concat, self, otherMat) )

    def qr( self, ROnly=0 ):
        'QR decomposition using Householder reflections: Q*R==self, Q.tr()*Q==I(n), R upper triangular'
        R = self
        m, n = R.size
        for i in range(min(m,n)):
            v, beta = R.tr()[i].house(i)
            R -= v.outer( R.tr().mmul(v)*beta )
        for i in range(1,min(n,m)): R[i][:i] = [0] * i
        R = new_matrix(R[:n])
        if ROnly:
            return R
        Q = R.tr().solve(self.tr()).tr()       # Rt Qt = At    nn  nm  = nm
        self.qr = lambda r=0, c=`self`: not r and c==`self` and (Q,R) or Matrix.qr(self,r) #Cache result
        assert NPOST or m>=n and Q.size==(m,n) and isinstance(R,UpperTri) or m<n and Q.size==(m,m) and R.size==(m,n)
        assert NPOST or Q.mmul(R)==self and Q.tr().mmul(Q)==identity_matrix(min(m,n))
        return Q, R

    def _solve( self, b ):
        '''General matrices (incuding) are solved using the QR composition.
        For inconsistent cases, returns the least squares solution'''
        Q, R = self.qr()
        return R.solve( Q.tr().mmul(b) )

    def solve( self, b ):
        'Divide matrix into a column vector or matrix and iterate to improve the solution'
        if b.dim==2:
            return new_matrix( map(self.solve, b.tr()) ).tr()
        assert NPRE or self.rows == len(b), 'Matrix row count %d must match vector length %d' % (self.rows, len(b))
        x = self._solve( b )
        diff = b - self.mmul(x)
        maxdiff = diff.dot(diff)
        for i in range(10):
            xnew = x + self._solve( diff )
            diffnew = b - self.mmul(xnew)
            maxdiffnew = diffnew.dot(diffnew)
            if maxdiffnew >= maxdiff:  break
            x, diff, maxdiff = xnew, diffnew, maxdiffnew
            #print >> sys.stderr, i+1, maxdiff
        assert NPOST or self.rows!=self.cols or self.mmul(x) == b
        return x

    def rank( self ):
        return Vector([ not row.forall(iszero) for row in self.qr(ROnly=1) ]).sum()

    def row_means(self):
        """
        Return vector consisting of row means.
        """
        return Vector([float(sum(v))/len(v) for v in self])

    def center_rows(self):
        """
        Return new matrix consisting of each the difference of each cell from
        its row mean.
        """
        row_means = self.row_means()
        m = []
        for i, v in enumerate(self):
            m.append(v - row_means[i])
        return Matrix(m)

    def covariance_by_rows(self, population_variance=False):
        """
        Returns covariance matrix (variables by rows).
        """
        x = self.center_rows()
        if population_variance:
            n = len(x[0])
        else:
            n = len(x[0])-1
        return x.mmul(x.tr()) / len(x[0])

    def col_means(self):
        """
        Return vector consisting of column means.
        """
        return self.tr().row_means()

    def center_cols(self):
        """
        Return new matrix consisting of each the difference of each cell from
        its column mean.
        """
        return self.tr().center_rows().tr()

    def covariance_by_cols(self, population_variance=False):
        """
        Returns covariance matrix (variables by columns).
        """
        x = self.center_cols()
        if population_variance:
            n = len(x)
        else:
            n = len(x)-1
        return x.tr().mmul(x) / n

class Square(Matrix):

    def lu( self ):
        'Factor a square matrix into lower and upper triangular form such that L.mmul(U)==A'
        n = self.rows
        L, U = identity_matrix(n), new_matrix(self[:])
        for i in range(n):
            for j in range(i+1,U.rows):
                assert U[i][i] != 0.0, 'LU requires non-zero elements on the diagonal'
                L[j][i] = m = 1.0 * U[j][i] / U[i][i]
                U[j] -= U[i] * m
        assert NPOST or isinstance(L,LowerTri) and isinstance(U,UpperTri) and L*U==self
        return L, U

    def __pow__( self, exp ):
        'Raise a square matrix to an integer power (i.e. A**3 is the same as A.mmul(A.mmul(A))'
        assert NPRE or exp==int(exp) and exp>0, 'Matrix powers only defined for positive integers not %s' % exp
        if exp == 1:
            return self
        if exp&1:
            return self.mmul(self ** (exp-1))
        sqrme = self ** (exp/2)
        return sqrme.mmul(sqrme)

    def det( self ):
        return self.qr( ROnly=1 ).det()

    def inverse( self ):
        return self.solve( identity_matrix(self.rows) )

    def hessenberg( self ):
        '''Householder reduction to Hessenberg Form (zeroes below the diagonal)
        while keeping the same eigenvalues as self.'''
        for i in range(self.cols-2):
            v, beta = self.tr()[i].house(i+1)
            self -= v.outer( self.tr().mmul(v)*beta )
            self -= self.mmul(v).outer(v*beta)
        return self

    def eigs( self ):
        'Estimate principal eigenvalues using the QR with shifts method'
        origTrace, origDet = self.trace(), self.det()
        self = self.hessenberg()
        eigvals = Vector([])
        for i in range(self.rows-1,0,-1):
            while not self[i][:i].forall(iszero):
                shift = identity_matrix(i+1) * self[i][i]
                q, r = (self - shift).qr()
                self = r.mmul(q) + shift
            eigvals.append( self[i][i] )
            self = new_matrix( [self[r][:i] for r in range(i)] )
        eigvals.append( self[0][0] )
        assert NPOST or iszero( (abs(origDet) - abs(eigvals.prod())) / 1000.0 )
        assert NPOST or iszero( origTrace - eigvals.sum() )
        return Vector(eigvals)

class Triangular(Square):

    def eigs( self ):
        return self.diag()

    def det( self ):
        return self.diag().prod()

class UpperTri(Triangular):

    def _solve( self, b ):
        'Solve an upper triangular matrix using backward substitution'
        x = Vector([])
        for i in range(self.rows-1, -1, -1):
            assert NPRE or self[i][i], 'Backsub requires non-zero elements on the diagonal'
            x.insert(0, (b[i] - x.dot(self[i][i+1:])) / self[i][i] )
        return x

class LowerTri(Triangular):

    def _solve( self, b ):
        'Solve a lower triangular matrix using forward substitution'
        x = Vector([])
        for i in range(self.rows):
            assert NPRE or self[i][i], 'Forward sub requires non-zero elements on the diagonal'
            x.append( (b[i] - x.dot(self[i][:i])) / self[i][i] )
        return x

def new_matrix( elems ):
    'Factory function to create a new matrix.'
    m, n = len(elems), len(elems[0])
    if m != n:
        return Matrix(elems)
    if n <= 1:
        return Square(elems)
    for i in range(1, len(elems)):
        if not iszero( max(map(abs, elems[i][:i])) ):
            break
    else:
        return UpperTri(elems)
    for i in range(0, len(elems)-1):
        if not iszero( max(map(abs, elems[i][i+1:])) ):
            return Square(elems)
    return LowerTri(elems)

def funToVector( tgtfun, low=-1, high=1, steps=40, EqualSpacing=0 ):
    '''Compute x,y points from evaluating a target function over an interval
    (low to high) at evenly spaces points or with Chebyshev abscissa spacing (
    default) '''
    if EqualSpacing:
        h = (0.0+high-low)/steps
        xvec = [low+h/2.0+h*i for i in range(steps)]
    else:
        scale, base = (0.0+high-low)/2.0, (0.0+high+low)/2.0
        xvec = [base+scale*math.cos(((2*steps-1-2*i)*math.pi)/(2*steps)) for i in range(steps)]
    yvec = map(tgtfun, xvec)
    return new_matrix( [xvec, yvec] )

def funfit( (xvec, yvec), basisfuns ):
    'Solves design matrix for approximating to basis functions'
    return new_matrix([ map(form,xvec) for form in basisfuns ]).tr().solve(Vector(yvec))

def polyfit( (xvec, yvec), degree=2 ):
    'Solves Vandermonde design matrix for approximating polynomial coefficients'
    return new_matrix([ [x**n for n in range(degree,-1,-1)] for x in xvec ]).solve(Vector(yvec))

def ratfit( (xvec, yvec), degree=2 ):
    'Solves design matrix for approximating rational polynomial coefficients (a*x**2 + b*x + c)/(d*x**2 + e*x + 1)'
    return new_matrix([[x**n for n in range(degree,-1,-1)]+[-y*x**n for n in range(degree,0,-1)] for x,y in zip(xvec,yvec)]).solve(Vector(yvec))

def genmat(m, n, func):
    if not n: n=m
    return new_matrix([ [func(i,j) for i in range(n)] for j in range(m) ])

def zeroes(m=1, n=None):
    'Zero matrix with side length m-by-m or m-by-n.'
    return genmat(m,n, lambda i,j: 0)

def identity_matrix(m=1, n=None):
    'Identity matrix with side length m-by-m or m-by-n'
    return genmat(m,n, lambda i,j: i==j)

def hilb(m=1, n=None):
    'Hilbert matrix with side length m-by-m or m-by-n.  Elem[i][j]=1/(i+j+1)'
    return genmat(m,n, lambda i,j: 1.0/(i+j+1.0))

def rand(m=1, n=None):
    'Random matrix with side length m-by-m or m-by-n'
    return genmat(m,n, lambda i,j: random.random())

if __name__ == '__main__':
    import cmath
    a = Table([1+2j,2,3,4])
    b = Table([5,6,7,8])
    C = Table([a,b])
    print 'a+b', a+b
    print '2+a', 2+a
    print 'a/5.0', a/5.0
    print '2*a+3*b', 2*a+3*b
    print 'a+C', a+C
    print '3+C', 3+C
    print 'C+b', C+b
    print 'C.sum()', C.sum()
    print 'C.map(math.cos)', C.map(cmath.cos)
    print 'C.conjugate()', C.conjugate()
    print 'C.real()', C.real()

    print zeroes(3)
    print identity_matrix(4)
    print hilb(3,5)

    C = new_matrix( [[1,2,3], [4,5,1,], [7,8,9]] )
    print C.mmul( C.tr()), '\n'
    print C ** 5, '\n'
    print C + C.tr(), '\n'

    A = C.tr().augment( new_matrix([[10,11,13]]).tr() ).tr()
    q, r = A.qr()
    assert q.mmul(r) == A
    assert q.tr().mmul(q)==identity_matrix(3)
    print 'q:\n', q, '\nr:\n', r, '\nQ.tr()&Q:\n', q.tr().mmul(q), '\nQ*R\n', q.mmul(r), '\n'
    b = Vector([50, 100, 220, 321])
    x = A.solve(b)
    print 'x:  ', x
    print 'b:  ', b
    print 'Ax: ', A.mmul(x)

    inv = C.inverse()
    print '\ninverse C:\n', inv, '\nC * inv(C):\n', C.mmul(inv)
    assert C.mmul(inv) == identity_matrix(3)

    points = (xvec,yvec) = funToVector(lambda x: math.sin(x)+2*math.cos(.7*x+.1), low=0, high=3, EqualSpacing=1)
    basis = [lambda x: math.sin(x), lambda x: math.exp(x), lambda x: x**2]
    print 'Func coeffs:', funfit( points, basis )
    print 'Poly coeffs:', polyfit( points, degree=5 )
    points = (xvec,yvec) = funToVector(lambda x: math.sin(x)+2*math.cos(.7*x+.1), low=0, high=3)
    print 'Rational coeffs:', ratfit( points )

    print polyfit(([1,2,3,4], [1,4,9,16]), 2)

    mtable = Vector([1,2,3]).outer(Vector([1,2]))
    print mtable, mtable.size

    A = new_matrix([ [2,0,3], [1,5,1], [18,0,6] ])
    print 'A:'
    print A
    print 'eigs:'
    print A.eigs()
    print 'Should be:', Vector([11.6158, 5.0000, -3.6158])
    print 'det(A)'
    print A.det()

    c = new_matrix( [[1,2,30],[4,5,10],[10,80,9]] )     # Failed example from Konrad Hinsen
    print 'C:\n', c
    print c.eigs()
    print 'Should be:', Vector([-8.9554, 43.2497, -19.2943])

    A = new_matrix([ [1,2,3,4], [4,5,6,7], [2,1,5,0], [4,2,1,0] ] )    # Kincaid and Cheney p.326
    print 'A:\n', A
    print A.eigs()
    print 'Should be:', Vector([3.5736, 0.1765, 11.1055, -3.8556])

    A = rand(3)
    q,r = A.qr()
    s,t = A.qr()
    print q is s                # Test caching
    print r is t
    A[1][1] = 1.1               # Invalidate the cache
    u,v = A.qr()
    print q is u                # Verify old result not used
    print r is v
    print u.mmul(v) == A        # Verify new result

    print 'Test qr on 3x5 matrix'
    a = rand(3,5)
    q,r = a.qr()
    print q.mmul(r) == a
    print q.tr().mmul(q) == identity_matrix(3)


