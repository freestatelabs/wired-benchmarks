""" Wired.py 
    Python implementation of the `Wired.jl` code, for benchmarking. 

    Optimization levels:
    0 :: Pure Python, standard library only (no vectorization or NumPy)
    1 :: NumPy
    2 :: Vectorization with NumPy 
    3 :: Vectorization with NumPy + Numba
    4 :: Parallel processing (on local CPU) via Numba
"""

from typing import Union
from math import sqrt, pi
import numpy as np 
from numba import njit, config, get_num_threads, set_num_threads, prange

Number = Union[int, float, np.number]
norm2 = np.linalg.norm

mu0 = 4*pi*(1e-7)

#
# Create the problem 
#

def create_problem(N: int, M: int, I=5000, Px=3, Py=3): 
    """ Make the sample problem: semi-infinite wire along the z-axis carrying 
        5 kA, output nodes are at x=y=3, along the z-axis. 

        N is the number of output nodes, M is the number input sources. Returns
        (nodes, sources), where nodes and sources are lists of lists.
    """

    # Make nodes
    zstart = 0; zend = 1
    dz = (zend - zstart) / (N - 1)
    z = zstart

    nodes = []
    for i in range(0, N):
        nodes.append([Px,Py,z])
        z += dz 

    # Make sources
    zstart = -1000; zend = 1000 
    dz = (zend - zstart) / N 
    z = zstart

    x = y = 0
    sources = [] 
    for i in range(0, M):
        sources.append([x, y, z, x, y, z+dz, I])
        z += dz 

    return nodes, sources


def create_problem2(N: int, M: int, I=5000, Px=3, Py=3): 
    """ Numpy version of the problem. """

    nodes, sources = create_problem(N, M, I, Px, Py)
    return np.array(nodes), np.array(sources)

def create_problem_rm(N: int, M: int, I=5000, Px=3, Py=3): 
    """ Numpy version of the problem, row-major. """

    nodes, sources = create_problem(N, M, I, Px, Py)
    return np.transpose(np.array(nodes)), np.array(sources)


#
# Optimization level 0 - pure Python 
# 

def cross(a: list, b: list):
    """ Cross-product for non-numpy implementation. """

    c1 = a[1]*b[2] - a[2]*b[1]
    c2 = a[2]*b[0] - a[0]*b[2] 
    c3 = a[0]*b[1] - a[1]*b[0]

    return [c1, c2, c3]


def norm(a) -> Number: 
    """ Vector norm for non-numpy implementation. """

    sumsq = 0 
    for a_ in a:
        sumsq += a_**2 

    return sqrt(sumsq)


def dot(a, b) -> Number: 
    """ Dot product of two vectors, for non-numpy implementation. """
    
    c = 0 
    for i in range(0, len(a)):
        c += a[i]*b[i] 

    return c


def minus(a: list, b:list) -> list:
    """ Subtract two lists from each other, for non-numpy implementation. """

    c = []
    for i in range(0, len(a)): 
        c.append(a[i] - b[i])

    return c


def biotsavart(J: list, K: list, I: Number, P: list) -> list:
    """ Opt level 0: compute the magnetic field caused by a single source at a
        single node. 
    """
    a = minus(K, J) 
    b = minus(J, P)
    c = minus(K, P)

    cxa = cross(c,a) 
    s0 = mu0*I/(4*pi) 
    s1 = norm(cxa)**2
    s2 = dot(a,c) 
    s3 = norm(c)
    s4 = dot(a,b)
    s5 = norm(b)

    bfield = cxa 
    for i in range(0, len(bfield)):
        bfield[i] *= s0 * (1 / s1) * (s2 / s3 - s4/s5)

    return bfield


def solve0(nodes: list, sources: list) -> list: 
    """ Optimization level 0: pure Python
    """

    Bfield = [] 

    for node in nodes: 

        b = [0, 0, 0]

        for source in sources: 
            b_ = biotsavart(source[0:3], source[3:6], source[6], node)

            for i in range(0, len(b)):
                b[i] += b_[i]

        Bfield.append(b)

    return Bfield


# 
# Optimization level 1: NumPy 
# 

def biotsavart1(node: np.ndarray, source: np.ndarray) -> np.ndarray: 
    """ Opt level 1: compute magnetic field at a node due to a single source
    """

    a = source[3:6] - source[0:3] 
    b = source[0:3] - node
    c = source[3:6] - node

    cxa = np.cross(c,a) 
    return (mu0*source[6]/(4*pi)) * cxa * (norm2(cxa)**-2) * (np.dot(a,c)/norm2(c)- np.dot(a,b,)/norm(b))


def solve1(nodes: np.ndarray, sources: np.ndarray) -> np.ndarray:
    """ Opt-level 1
    """ 

    Bfield = np.zeros(nodes.shape) 

    for source in sources: 
        for i in range(0, nodes.shape[0]):
            Bfield[i,:] += biotsavart1(nodes[i,:], source)

    return Bfield 


#
# Optimization level 2: Numpy vectorization
#

def normrows(A: np.ndarray) -> np.ndarray: 
    """ opt-level 2: Take the norm of the rows of A 
        Required for numba in opt-level 3
    """
    return np.sqrt(A[:,0]**2 + A[:,1]**2 + A[:,2]**2)


def dotrows(A: np.ndarray, b: np.ndarray): 
    """c = A*b """
    return np.sum(A * b, axis=1)


def biotsavart2(nodes: np.ndarray, source: np.ndarray) -> np.ndarray: 
    """ Vectorized version
    """
    a = source[3:6] - source[0:3] 
    b = np.subtract(source[0:3], nodes)
    c = np.subtract(source[3:6], nodes)

    cxa = np.cross(c, a) 
    coeffs = (mu0*source[6]/(4*pi) / normrows(cxa)**2) * (dotrows(c,a) / normrows(c) - dotrows(b,a) / normrows(b))

    return cxa * coeffs[:,None]


def solve2(nodes, sources): 
    """ Vectorized using NumPy
    """

    Bfield = np.zeros(nodes.shape) 

    for source in sources: 
        Bfield += biotsavart2(nodes, source)

    return Bfield

#
# Optimization level 2: Numpy vectorization (row-major)
# Doesn't seem to be any faster or slower than the default 
#

def normcols(A: np.ndarray) -> np.ndarray: 
    """ opt-level 2: Take the norm of the rows of A 
        Required for numba in opt-level 3
    """
    return np.sqrt(A[0,:]**2 + A[1,:]**2 + A[2,:]**2)


def dotcols(A: np.ndarray, b: np.ndarray): 
    """c = A*b """
    return np.sum(A * np.reshape(b,(3,1)), axis=0)


def biotsavart2_rm(nodes: np.ndarray, source: np.ndarray) -> np.ndarray: 
    """ Vectorized version
    """
    a = source[3:6] - source[0:3] 
    b = -nodes + np.reshape(source[0:3],(3,1)) #np.subtract(np.reshape(source[0:3],(1,3)), nodes)
    c = -nodes + np.reshape(source[3:6],(3,1))

    cxa = np.reshape(np.cross(c, a, axisa=0),(3,nodes.shape[1]))
    coeffs = (mu0*source[6]/(4*pi) / normcols(cxa)**2) * (dotcols(c,a) / normcols(c) - dotcols(b,a) / normcols(b))

    return cxa * np.reshape(coeffs[:,None], (1, nodes.shape[1]))


def solve2_rm(nodes, sources): 
    """ Vectorized using NumPy
    """

    Bfield = np.zeros(nodes.shape) 

    for source in sources: 
        Bfield += biotsavart2_rm(nodes, source)

    return Bfield


# 
# Optimization level 3: Numba 
#  

nbnorm = njit(normrows)
nbdot = njit(dotrows)


@njit
def solve3(nodes: np.ndarray, sources: np.ndarray) -> np.ndarray: 

    Ns = sources.shape[0]
    Bfield = np.zeros(nodes.shape) 

    for i in range(0, Ns): 
        source = sources[i,:]
        a = source[3:6] - source[0:3] 
        b = np.subtract(source[0:3], nodes)
        c = np.subtract(source[3:6], nodes)

        cxa = np.cross(c, a) 
        cxa *= (mu0*source[6]/(4*pi)) * (1 / nbnorm(cxa)**2) * (nbdot(c,a) / nbnorm(c) - nbdot(b,a) / nbnorm(b))
        Bfield = np.add(Bfield, cxa)
    return Bfield


#
# Optimization level 4: parallel Numba
#

@njit
def threadindices(it, Nt, N):

        Nperthread, remainder = divmod(N, Nt)

        if it == 1:
            i1 = 1 
            i2 = i1 + Nperthread - 1 
        elif it == Nt:
            i2 = N 
            i1 = i2 - Nperthread - remainder + 1 
        else:
            i1 = (it-1)*Nperthread + 1
            i2 = i1 + Nperthread - 1 

        return i1,i2

# Íconfig.THREADING_LAYER = 'threadsafe'
@njit(parallel=True)
def solve4(nodes: np.ndarray, sources: np.ndarray):

    Bfield = np.zeros(nodes.shape) 
    Nt = get_num_threads()
    Ns = sources.shape[0]
    Bfield_ = np.zeros((Nt, nodes.shape[0], nodes.shape[1]) )

    for i in prange(0, Nt):

        i1, i2 = threadindices(i, Nt, Ns)

        for j in range(i1, i2+1): #[i1:i2,:]: 
            source = sources[j,:]
            a = source[3:6] - source[0:3] 
            b = np.subtract(source[0:3], nodes)
            c = np.subtract(source[3:6], nodes)

            cxa = np.cross(c, a) # Nx3 matrix
            f = (mu0*source[6]/(4*pi) / nbnorm(cxa)**2) * (nbdot(c,a) / nbnorm(c) - nbdot(b,a) / nbnorm(b))
            cxa[:,0] *= f 
            cxa[:,1] *= f 
            cxa[:,2] *= f
            Bfield_[i] = np.add(Bfield_[i], cxa) 
        
    Bfield = np.sum(Bfield_, axis=0)

    return Bfield

# For convenience:
solvers = [solve0, solve1, solve2, solve3, solve4]