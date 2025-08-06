# pyright: reportOperatorIssue=false
# pyright: reportArgumentType=false

from sympy.combinatorics.permutations import Permutation
import sympy as sy

# from sympy.matrices import sy.Matrix
# from sympy.functions import sy.sin, sy.cos
# from sympy.core.numbers import sy.I
# from sympy.functions.elementary.exponential import sy.exp
# from sympy.core import Expr
#from typing import List, Dict, Any, Literal


'''All book references in this file are related to
Professor Christoph Berger book 
"Particle Physics: From the Basics to Modern Experiments",
ISBN-13:	9783662717059'''

#
Eins: sy.Matrix = sy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
#
g0: sy.Matrix = sy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
g1: sy.Matrix = sy.Matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]])
g2: sy.Matrix = sy.Matrix([[0, 0, 0, -sy.I], [0, 0, sy.I, 0],
             [0, sy.I, 0, 0], [-sy.I, 0, 0, 0]])
g3: sy.Matrix = sy.Matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]])

# g5=sy.I*g0*g1*g2*g3
g5: sy.Matrix = sy.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
one: sy.Matrix = sy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


proj_r = (one + g5) / 2
proj_l = (one - g5) / 2
#
# Pauli Matrices (compare book chapter 2.2 or 2.8)
#
sig1 = sy.Matrix([[0, 1], [1, 0]])
sig2 = sy.Matrix([[0, -sy.I], [sy.I, 0]])
sig3 = sy.Matrix([[1, 0], [0, -1]])


def chi1(theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    
    elm1 = sy.cos(theta / 2) # type: ignore
    elm2 = sy.exp(sy.I * phi) * sy.sin(theta / 2)  # type: ignore
    result = sy.Matrix([[elm1], [elm2]])
    return result


def chi2(theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm1: sy.Expr = -sy.exp(-sy.I * phi) * sy.sin(theta / 2)  # type: ignore
    elm2: sy.Expr = sy.cos(theta / 2)  # type: ignore
    result = sy.Matrix([[elm1], [elm2]])
    return result


#
# Helicity-spinors (compare book chapter 3.1)
#
def u_min(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm1 = -sy.exp(-sy.I * phi) * sy.sin(theta / 2) # type: ignore
    elm2 = sy.cos(theta / 2) # type: ignore
    fac1 = sy.sqrt(Ef - mf) / sy.sqrt(Ef + mf) # type: ignore
    fac2 = sy.sqrt(Ef + mf) # type: ignore
    elm3 = -fac1 * elm1
    elm4 = -fac1 * elm2
    result = fac2 * sy.Matrix([[elm1], [elm2], [elm3], [elm4]])
    return result


#
def u_plus(Ef: sy.Symbol, mf: sy.Symbol, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm2 = sy.exp(sy.I * phi) * sy.sin(theta / 2) # type: ignore
    elm1 = sy.cos(theta / 2) # type: ignore
    fac1 = sy.sqrt(Ef - mf) / sy.sqrt(Ef + mf) # type: ignore
    fac2 = sy.sqrt(Ef + mf) # type: ignore
    elm3 = fac1 * elm1
    elm4 = fac1 * elm2
    result = fac2 * sy.Matrix([[elm1], [elm2], [elm3], [elm4]])
    return result


#
def v_min(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm4 = sy.exp(sy.I * phi) * sy.sin(theta / 2) # pyright: ignore[reportOperatorIssue]
    elm3 = sy.cos(theta / 2) 
    fac1 = sy.sqrt(Ef - mf) / sy.sqrt(Ef + mf) 
    fac2 = sy.sqrt(Ef + mf)
    elm1 = fac1 * elm3
    elm2 = fac1 * elm4
    result = fac2 * sy.Matrix([[elm1], [elm2], [elm3], [elm4]])
    return result


#
def v_plus(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm3 = -sy.exp(-sy.I * phi) * sy.sin(theta / 2) # type: ignore
    elm4 = sy.cos(theta / 2) # type: ignore
    fac1 = sy.sqrt(Ef - mf) / sy.sqrt(Ef + mf) # type: ignore
    fac2 = sy.sqrt(Ef + mf) # type: ignore
    elm1 = -fac1 * elm3
    elm2 = -fac1 * elm4
    result = fac2 * sy.Matrix([[elm1], [elm2], [elm3], [elm4]])
    return result


#
def u_minbar(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    h1 = u_min(Ef, mf, theta, phi).T * g0 # type: ignore
    elm1 = h1[0] * sy.exp(2 * sy.I * phi) # type: ignore
    elm2 = h1[1]
    # fac1=sy.sqrt(Ef-mf)/sy.sqrt(Ef+mf)
    # fac2=sy.sqrt(Ef+mf)
    elm3 = h1[2] * sy.exp(2 * sy.I * phi)
    elm4 = h1[3]
    result = sy.Matrix([[elm1, elm2, elm3, elm4]])
    return result


#
def u_plusbar(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    h1 = u_plus(Ef, mf, theta, phi).T * g0
    elm2 = h1[1] * sy.exp(-2 * sy.I * phi)
    elm1 = h1[0]
    # fac1=sy.sqrt(Ef-mf)/sy.sqrt(Ef+mf)
    # fac2=sy.sqrt(Ef+mf)
    elm4 = h1[3] * sy.exp(-2 * sy.I * phi)
    elm3 = h1[2]
    result = sy.Matrix([[elm1, elm2, elm3, elm4]])
    return result


#
def v_minbar(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    h1 = v_min(Ef, mf, theta, phi).T * g0
    elm2 = h1[1] * sy.exp(-2 * sy.I * phi)
    elm1 = h1[0]
    # fac1=sy.sqrt(Ef-mf)/sy.sqrt(Ef+mf)
    # fac2=sy.sqrt(Ef+mf)
    elm4 = h1[3] * sy.exp(-2 * sy.I * phi)
    elm3 = h1[2]
    result = sy.Matrix([[elm1, elm2, elm3, elm4]])
    return result


#
def v_plusbar(Ef: sy.Expr, mf: sy.Expr, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    h1 = v_plus(Ef, mf, theta, phi).T * g0
    elm1 = h1[0] * sy.exp(2 * sy.I * phi)
    elm2 = h1[1]
    # fac1=sy.sqrt(Ef-mf)/sy.sqrt(Ef+mf)
    # fac2=sy.sqrt(Ef+mf)
    elm3 = h1[2] * sy.exp(2 * sy.I * phi)
    elm4 = h1[3]
    result = sy.Matrix([[elm1, elm2, elm3, elm4]])
    return result


#
def u(pp: sy.Matrix, ha: int) -> sy.Matrix:
    '''This function defines the spinor $u(pe,1)$, that represents
    is a particle spinor with momentum $pp$ and positive helicity ($\lambda=1/2$), $v(pe,-1)$
    an antiparticle spinor with negative helicity ($\lambda=-1/2$).
    pp: fourmomenta of the particle
    ha: helicity (positive (1), null(0), or negative(-1))'''

    if ha >= 0:
        result = u_plus(pp[0], pp[1], pp[2], pp[3])
    elif ha < 0:
        result = u_min(pp[0], pp[1], pp[2], pp[3])
    return result


#
def v(pp, ha):
    '''This function defines the spinor $v(pe,1)$, that represents
    is a antiparticle spinor with momentum $pp$ and positive helicity ($\lambda=1/2$)
    pp: fourmomenta of the particle
    ha: helicity (positive (1), null(0), or negative(-1))'''

    if ha >= 0:
        result = v_plus(pp[0], pp[1], pp[2], pp[3])
    elif ha < 0:
        result = v_min(pp[0], pp[1], pp[2], pp[3])
    return result


#
def ubar(pp, ha):
    '''
    pp: fourmomenta of the particle
    ha: helicity (positive (1), null(0), or negative(-1))
    '''
    if ha >= 0:
        result = u_plusbar(pp[0], pp[1], pp[2], pp[3])
    elif ha < 0:
        result = u_minbar(pp[0], pp[1], pp[2], pp[3])
    return result


#
def vbar(pp, ha):
    '''
    pp: fourmomenta of the particle
    ha: helicity (positive (1), null(0), or negative(-1))
    '''
    if ha >= 0:
        result = v_plusbar(pp[0], pp[1], pp[2], pp[3])
    elif ha < 0:
        result = v_minbar(pp[0], pp[1], pp[2], pp[3])
    return result


#
# Polarization vectors for photons and W bosons (compare book chapter 2.4
# and chapter 3.2)
#
def pol_plus(theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm1 = 0
    h1 = sy.cos(theta / 2) ** 2 - sy.sin(theta / 2) ** 2 # type: ignore
    h2 = 2 * sy.sin(theta / 2) * sy.cos(theta / 2) # type: ignore
    elm2 = sy.sqrt(2) * sy.cos(phi) * h1 / 2 - sy.I * sy.sqrt(2) * sy.sin(phi) / 2 # type: ignore
    elm3 = sy.sqrt(2) * sy.sin(phi) * h1 / 2 + sy.I * sy.sqrt(2) * sy.cos(phi) / 2 # type: ignore
    elm4 = -sy.sqrt(2) * h2 / 2
    result = [elm1, elm2, elm3, elm4]
    return sy.Matrix(result)


#
#
def pol_min(theta, phi) -> sy.Matrix:
    elm1 = 0
    h1 = sy.cos(theta / 2) ** 2 - sy.sin(theta / 2) ** 2 # type: ignore
    h2 = 2 * sy.sin(theta / 2) * sy.cos(theta / 2) # type: ignore
    elm2 = sy.sqrt(2) * sy.cos(phi) * h1 / 2 + sy.I * sy.sqrt(2) * sy.sin(phi) / 2 # type: ignore
    elm3 = sy.sqrt(2) * sy.sin(phi) * h1 / 2 - sy.I * sy.sqrt(2) * sy.cos(phi) / 2 # type: ignore
    elm4 = -sy.sqrt(2) * h2 / 2
    result = [elm1, elm2, elm3, elm4]
    return sy.Matrix(result)


#
#
def pol_0(e0, m0, theta: sy.Expr, phi: sy.Expr) -> sy.Matrix:
    elm1 = sy.sqrt(e0**2 - m0**2) / m0
    elm2 = e0 * sy.sin(theta) * sy.cos(phi) / m0
    elm3 = e0 * sy.sin(theta) * sy.sin(phi) / m0
    elm4 = e0 * sy.cos(theta) / m0
    result = [elm1, elm2, elm3, elm4]
    return sy.Matrix(result)


#
def pol(pp, ha):
    e0 = pp[0]
    m0 = pp[1]
    theta = pp[2]
    phi = pp[3]
    if ha > 0:
        result = pol_plus(theta, phi)
    elif ha < 0:
        result = pol_min(theta, phi)
    else:
        result = pol_0(e0, m0, theta, phi)

    return result


#
def polbar(pp, ha):
    e0 = pp[0]
    m0 = pp[1]
    theta = pp[2]
    phi = pp[3]
    if ha > 0:
        result = pol_min(theta, phi)
        # display(f'+ pol {result}')
    elif ha < 0:
        result = pol_plus(theta, phi)
        # display(f'- pol {result}')
    else:
        result = pol_0(e0, m0, theta, phi)
        # display(f'0 pol {result}')
    return sy.transpose(result)  # Modified 2023-12-13, add transpose conjugation


#
# useful routines for calculations:
#
# fourvec calculates E, p_x,p_y,p_z from E, m, theta, phi
#
def fourvec(pp):
    elm1 = pp[0]
    elm2 = pp[1]
    elm3 = pp[2]
    elm4 = pp[3]

    h0 = elm3 / 2
    h1 = 2 * sy.sin(h0) * sy.cos(h0) # type: ignore
    h2 = (sy.cos(h0) ** 2) - (sy.sin(h0) ** 2) # type: ignore
    p = sy.sqrt(elm1**2 - elm2**2)
    # result = [elm1, p * h1 * sy.cos(elm4), p * h1 * sy.sin(elm4), p * h2]
    result = sy.Matrix([elm1, p * h1 * sy.cos(elm4), p * h1 * sy.sin(elm4), p * h2])
    return result


#
#
# dotprod4 evaluates the scalar product of 2 fourvectors
#
def dotprod4(jj1, jj2):
    h1 = jj1[0] * jj2[0] - jj1[1] * jj2[1] - jj1[2] * jj2[2] - jj1[3] * jj2[3]
    result = sy.simplify(h1)
    return result


#
# dag evaluates the dagger matrix as defined in chapter 3.1
# NOTE: This is in fact what we call the "slash"
def dag(pp):
    elm1 = pp[0]
    elm2 = pp[1]
    elm3 = pp[2]
    elm4 = pp[3]
    result = elm1 * g0 - elm2 * g1 - elm3 * g2 - elm4 * g3
    return result


#
# currents gamma_mu coupling (QED, compare book chapter 3.2)
#
def ubv(p1, h1, p2, h2):
    j0 = ubar(p1, h1) * g0 * v(p2, h2)
    j1 = ubar(p1, h1) * g1 * v(p2, h2)
    j2 = ubar(p1, h1) * g2 * v(p2, h2)
    j3 = ubar(p1, h1) * g3 * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
def ubu(p1, h1, p2, h2):
    j0 = ubar(p1, h1) * g0 * u(p2, h2)
    j1 = ubar(p1, h1) * g1 * u(p2, h2)
    j2 = ubar(p1, h1) * g2 * u(p2, h2)
    j3 = ubar(p1, h1) * g3 * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
def vbu(p1, h1, p2, h2):
    j0 = vbar(p1, h1) * g0 * u(p2, h2)
    j1 = vbar(p1, h1) * g1 * u(p2, h2)
    j2 = vbar(p1, h1) * g2 * u(p2, h2)
    j3 = vbar(p1, h1) * g3 * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
def vbv(p1, h1, p2, h2):
    j0 = vbar(p1, h1) * g0 * v(p2, h2)
    j1 = vbar(p1, h1) * g1 * v(p2, h2)
    j2 = vbar(p1, h1) * g2 * v(p2, h2)
    j3 = vbar(p1, h1) * g3 * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
# charged current weak interaction, (1-g5) coupling (compare book chapter 6.1)
#
#
def ubvw(p1, h1, p2, h2):
    j0 = ubar(p1, h1) * g0 * proj_l * v(p2, h2)
    j1 = ubar(p1, h1) * g1 * proj_l * v(p2, h2)
    j2 = ubar(p1, h1) * g2 * proj_l * v(p2, h2)
    j3 = ubar(p1, h1) * g3 * proj_l * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
#
def ubuw(p1, h1, p2, h2):
    j0 = ubar(p1, h1) * g0 * proj_l * u(p2, h2)
    j1 = ubar(p1, h1) * g1 * proj_l * u(p2, h2)
    j2 = ubar(p1, h1) * g2 * proj_l * u(p2, h2)
    j3 = ubar(p1, h1) * g3 * proj_l * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
#
def vbuw(p1, h1, p2, h2):
    j0 = vbar(p1, h1) * g0 * proj_l * u(p2, h2)
    j1 = vbar(p1, h1) * g1 * proj_l * u(p2, h2)
    j2 = vbar(p1, h1) * g2 * proj_l * u(p2, h2)
    j3 = vbar(p1, h1) * g3 * proj_l * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
def vbvw(p1, h1, p2, h2):
    j0 = vbar(p1, h1) * g0 * proj_l * v(p2, h2)
    j1 = vbar(p1, h1) * g1 * proj_l * v(p2, h2)
    j2 = vbar(p1, h1) * g2 * proj_l * v(p2, h2)
    j3 = vbar(p1, h1) * g3 * proj_l * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
# currents cv-ca*g5 coupling (compare book chapter 6.4)
#
def ubvva(p1, h1, p2, h2, cv, ca):
    gv = cv / 2
    ga = ca / 2
    j0 = ubar(p1, h1) * g0 * (gv * one - ga * g5) * v(p2, h2)
    j1 = ubar(p1, h1) * g1 * (gv * one - ga * g5) * v(p2, h2)
    j2 = ubar(p1, h1) * g2 * (gv * one - ga * g5) * v(p2, h2)
    j3 = ubar(p1, h1) * g3 * (gv * one - ga * g5) * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
#
def ubuva(p1, h1, p2, h2, cv, ca):
    gv = cv / 2
    ga = ca / 2
    j0 = ubar(p1, h1) * g0 * (gv * one - ga * g5) * u(p2, h2)
    j1 = ubar(p1, h1) * g1 * (gv * one - ga * g5) * u(p2, h2)
    j2 = ubar(p1, h1) * g2 * (gv * one - ga * g5) * u(p2, h2)
    j3 = ubar(p1, h1) * g3 * (gv * one - ga * g5) * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


##
def vbuva(p1, h1, p2, h2, cv, ca):
    gv = cv / 2
    ga = ca / 2
    j0 = vbar(p1, h1) * g0 * (gv * one - ga * g5) * u(p2, h2)
    j1 = vbar(p1, h1) * g1 * (gv * one - ga * g5) * u(p2, h2)
    j2 = vbar(p1, h1) * g2 * (gv * one - ga * g5) * u(p2, h2)
    j3 = vbar(p1, h1) * g3 * (gv * one - ga * g5) * u(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
#
def vbvva(p1, h1, p2, h2, cv, ca):
    gv = cv / 2
    ga = ca / 2
    j0 = vbar(p1, h1) * g0 * (gv * one - ga * g5) * v(p2, h2)
    j1 = vbar(p1, h1) * g1 * (gv * one - ga * g5) * v(p2, h2)
    j2 = vbar(p1, h1) * g2 * (gv * one - ga * g5) * v(p2, h2)
    j3 = vbar(p1, h1) * g3 * (gv * one - ga * g5) * v(p2, h2)
    result = [sy.simplify(j0[0]), sy.simplify(
        j1[0]), sy.simplify(j2[0]), sy.simplify(j3[0])]
    return result


#
#


#
# Pauli currents as needed in book chapter 5.2
#
def ubuP(p1, h1, p2, h2):
    M = p1[1] # type: ignore
    p14 = fourvec(p1) 
    p24 = fourvec(p2)
    # qmu with lower index
    q0 = p14[0] - p24[0] # type: ignore
    q1 = -p14[1] + p24[1] # type: ignore
    q2 = -p14[2] + p24[2] # type: ignore
    q3 = -p14[3] + p24[3] # type: ignore
    # sigmamunu with upper index like in Bjorken Drell
    # but without factor i/2
    jF0 = (g0 * g1 - g1 * g0) * q1 + (g0 * g2 - g2 * g0) * \
        q2 + (g0 * g3 - g3 * g0) * q3
    jF1 = (g1 * g0 - g0 * g1) * q0 + (g1 * g2 - g2 * g1) * \
        q2 + (g1 * g3 - g3 * g1) * q3
    jF2 = (g2 * g0 - g0 * g2) * q0 + (g2 * g1 - g1 * g2) * \
        q1 + (g2 * g3 - g3 * g2) * q3
    jF3 = (g3 * g0 - g0 * g3) * q0 + (g3 * g1 - g1 * g3) * \
        q1 + (g3 * g2 - g2 * g3) * q2
    j0 = ubar(p1, h1) * jF0 / 4 / M * u(p2, h2)
    j1 = ubar(p1, h1) * jF1 / 4 / M * u(p2, h2)
    j2 = ubar(p1, h1) * jF2 / 4 / M * u(p2, h2)
    j3 = ubar(p1, h1) * jF3 / 4 / M * u(p2, h2)
    # include factor i/2
    result = [-j0[0], -j1[0], -j2[0], -j3[0]]
    return result

#
# s channel and u channel amplitudes for compton scattering
#
def compts(k1, hk1, p1, hp1, k2, hk2, p2, hp2):
    ke = fourvec(k1)
    pe = fourvec(p1)
    mm = p1[1]
    mk = k1[1]
    nenner = 2 * dotprod4(ke, pe) + mk**2
    #        nenner=s0+mk^2
    epsbar = polbar(k2, hk2)
    eps = pol(k1, hk1)
    kern = dag(epsbar) * (dag(ke) + dag(pe) + mm * one) * dag(eps)
    h1 = sy.simplify(ubar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0] / nenner)
    return result


#
def comptu(k1, hk1, p1, hp1, k2, hk2, p2, hp2):
    ka = fourvec(k2)
    pe = fourvec(p1)
    mm = p1[1]
    mk = k2[1]
    nenner = -2 * dotprod4(ka, pe) + mk**2
    # nenner:=u0+mk^2:
    epsbar = polbar(k2, hk2)
    eps = pol(k1, hk1)
    kern = dag(eps) * (dag(pe) - dag(ka) + mm * one) * dag(epsbar)
    h1 = sy.simplify(ubar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0] / nenner)
    return result


#
def compt(k1, hk1, p1, hp1, k2, hk2, p2, hp2):
    t1 = compts(k1, hk1, p1, hp1, k2, hk2, p2, hp2)
    t2 = comptu(k1, hk1, p1, hp1, k2, hk2, p2, hp2)
    result = sy.simplify(t1 + t2)
    return result


#
#
# Sometimes one wants the amplitudes without the propagator term
#
def Ncompts(k1, hk1, p1, hp1, k2, hk2, p2, hp2):
    ke = fourvec(k1)
    pe = fourvec(p1)
    mm = p1[1]
    mk = k1[1]
    epsbar = polbar(k2, hk2)
    eps = pol(k1, hk1)
    kern = dag(epsbar) * (dag(ke) + dag(pe) + mm * one) * dag(eps)
    h1 = sy.simplify(ubar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0])
    return result


#
def Ncomptu(k1, hk1, p1, hp1, k2, hk2, p2, hp2):
    ka = fourvec(k2)
    pe = fourvec(p1)
    mm = p1[1]
    mk = k1[1]
    epsbar = polbar(k2, hk2)
    eps = pol(k1, hk1)
    kern = dag(eps) * (dag(pe) - dag(ka) + mm * one) * dag(epsbar)
    h1 = sy.simplify(ubar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0])
    return result


#
#
# Now f+fbar-->gamma gamma
# p1+p2=k1+k2
#
#
def gamgamt(p1, hp1, p2, hp2, k1, hk1, k2, hk2):
    ka = fourvec(k1)
    pe = fourvec(p1)
    mm = p1[1]
    M = k1[1]
    nenner = -2 * dotprod4(ka, pe) + M**2
    epsbar2 = polbar(k2, hk2)
    epsbar1 = polbar(k1, hk1)
    kern = dag(epsbar2) * (dag(ka) - dag(pe) + mm * one) * dag(epsbar1)
    h1 = sy.simplify(vbar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0] / nenner)
    return result


#
#
def gamgamu(p1, hp1, p2, hp2, k1, hk1, k2, hk2):
    ka = fourvec(k2)
    pe = fourvec(p1)
    mm = p1[1]
    M = k2[1]
    nenner = -2 * dotprod4(ka, pe) + M**2
    epsbar2 = polbar(k2, hk2)
    epsbar1 = polbar(k1, hk1)
    kern = dag(epsbar1) * (dag(ka) - dag(pe) + mm * one) * dag(epsbar2)
    h1 = sy.simplify(vbar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0] / nenner)
    return result


#
#
def gamgam(p1, hp1, p2, hp2, k1, hk1, k2, hk2):
    t1 = gamgamt(p1, hp1, p2, hp2, k1, hk1, k2, hk2)
    t2 = gamgamu(p1, hp1, p2, hp2, k1, hk1, k2, hk2)
    result = sy.simplify(t1 + t2)
    return result


#
# Amplitudes without denominator
#
#
def Ngamgamt(p1, hp1, p2, hp2, k1, hk1, k2, hk2):
    ka = fourvec(k1)
    pe = fourvec(p1)
    mm = p1[1]
    epsbar2 = polbar(k2, hk2)
    epsbar1 = polbar(k1, hk1)
    kern = dag(epsbar2) * (dag(ka) - dag(pe) + mm * one) * dag(epsbar1)
    h1 = sy.simplify(vbar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0])
    return result


#
#
def Ngamgamu(p1, hp1, p2, hp2, k1, hk1, k2, hk2):
    ka = fourvec(k2)
    pe = fourvec(p1)
    mm = p1[1]
    epsbar2 = polbar(k2, hk2)
    epsbar1 = polbar(k1, hk1)
    kern = dag(epsbar1) * (dag(ka) - dag(pe) + mm * one) * dag(epsbar2)
    h1 = sy.simplify(vbar(p2, hp2) * kern * u(p1, hp1))
    result = sy.simplify(h1[0, 0])
    return result


##
#
# The routines V3g... and V4g... are needed for calculating the cross section
# of elastic gluon gluon scattering, gg to gg as treated in chapter 4.2 of the book
# and in the noteboo ggtogg.ipynb
#
#
def V3gtchannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    r1 = V3gInOut(p1, hp1, p3, hp3)
    r2 = V3gInOut(p2, hp2, p4, hp4)
    result = dotprod4(r1, r2)
    return result


#
#
def V4gtchannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    eps1 = pol(p1, hp1)
    eps2 = pol(p2, hp2)
    eps3 = polbar(p3, hp3)
    eps4 = polbar(p4, hp4)
    r1 = dotprod4(eps1, eps4)
    r2 = dotprod4(eps2, eps3)
    r3 = dotprod4(eps1, eps2)
    r4 = dotprod4(eps3, eps4)
    result = -r1 * r2 + r3 * r4
    return result


#
def V3guchannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    r1 = V3gInOut(p1, hp1, p4, hp4)
    r2 = V3gInOut(p2, hp2, p3, hp3)
    result = dotprod4(r1, r2)
    return result


#
def V4guchannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    eps1 = pol(p1, hp1)
    eps2 = pol(p2, hp2)
    eps3 = polbar(p3, hp3)
    eps4 = polbar(p4, hp4)
    r1 = dotprod4(eps1, eps2)
    r2 = dotprod4(eps3, eps4)
    r3 = dotprod4(eps1, eps3)
    r4 = dotprod4(eps2, eps4)
    result = r1 * r2 - r3 * r4
    return result


#
#
def V3gschannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    r1 = V3gInIn(p1, hp1, p2, hp2)
    r2 = V3gOutOut(p3, hp3, p4, hp4)
    result = dotprod4(r1, r2)
    return result


#
def V4gschannel(p1, hp1, p2, hp2, p3, hp3, p4, hp4):
    eps1 = pol(p1, hp1)
    eps2 = pol(p2, hp2)
    eps3 = polbar(p3, hp3)
    eps4 = polbar(p4, hp4)
    r1 = dotprod4(eps1, eps3)
    r2 = dotprod4(eps2, eps4)
    r3 = dotprod4(eps1, eps4)
    r4 = dotprod4(eps2, eps3)
    result = r1 * r2 - r3 * r4
    return result


#
#
# The  following routines, e.e. V3gInIn re useful for certain checks in undertsanding
# the notebook ggtogg.ipynb
#
def V3gInIn(k, hk, p, hp):
    eps1 = pol(k, hk)
    eps2 = pol(p, hp)
    p1 = fourvec(k)
    p2 = fourvec(p)
    s1 = dotprod4(eps1, eps2)
    s2 = dotprod4(eps1, p2)
    s3 = dotprod4(eps2, p1)
    elm1 = s1 * (p1[0] - p2[0]) + 2 * s2 * eps2[0] - 2 * s3 * eps1[0]
    elm2 = s1 * (p1[1] - p2[1]) + 2 * s2 * eps2[1] - 2 * s3 * eps1[1]
    elm3 = s1 * (p1[2] - p2[2]) + 2 * s2 * eps2[2] - 2 * s3 * eps1[2]
    elm4 = s1 * (p1[3] - p2[3]) + 2 * s2 * eps2[3] - 2 * s3 * eps1[3]
    result = [elm1, elm2, elm3, elm4]
    return result


def V3gOutOut(k, hk, p, hp):
    eps3 = polbar(k, hk)
    eps4 = polbar(p, hp)
    p3 = fourvec(k)
    p4 = fourvec(p)
    s1 = dotprod4(eps3, eps4)
    s2 = dotprod4(eps3, p4)
    s3 = dotprod4(eps4, p3)
    elm1 = s1 * (p4[0] - p3[0]) - 2 * s2 * eps4[0] + 2 * s3 * eps3[0]
    elm2 = s1 * (p4[1] - p3[1]) - 2 * s2 * eps4[1] + 2 * s3 * eps3[1]
    elm3 = s1 * (p4[2] - p3[2]) - 2 * s2 * eps4[2] + 2 * s3 * eps3[2]
    elm4 = s1 * (p4[3] - p3[3]) - 2 * s2 * eps4[3] + 2 * s3 * eps3[3]
    result = [elm1, elm2, elm3, elm4]
    return result


def V3gInOut(k, hk, p, hp):
    eps1 = pol(k, hk)
    eps3 = polbar(p, hp)
    p1 = fourvec(k)
    p3 = fourvec(p)
    s1 = dotprod4(eps1, eps3)
    s2 = dotprod4(eps1, p3)
    s3 = dotprod4(eps3, p1)
    elm1 = sy.simplify(s1 * (p1[0] + p3[0]) - 2 * s2 * eps3[0] - 2 * s3 * eps1[0])
    elm2 = sy.simplify(s1 * (p1[1] + p3[1]) - 2 * s2 * eps3[1] - 2 * s3 * eps1[1])
    elm3 = sy.simplify(s1 * (p1[2] + p3[2]) - 2 * s2 * eps3[2] - 2 * s3 * eps1[2])
    elm4 = sy.simplify(s1 * (p1[3] + p3[3]) - 2 * s2 * eps3[3] - 2 * s3 * eps1[3])
    result = [elm1, elm2, elm3, elm4]
    return result


#
# And of cours we need the SU3 Generators
#
lam1 = sy.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
lam2 = sy.Matrix([[0, -sy.I, 0], [sy.I, 0, 0], [0, 0, 0]])
lam3 = sy.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
lam4 = sy.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
lam5 = sy.Matrix([[0, 0, -sy.I], [0, 0, 0], [sy.I, 0, 0]])
lam6 = sy.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
lam7 = sy.Matrix([[0, 0, 0], [0, 0, -sy.I], [0, sy.I, 0]])
lam08 = sy.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]])
lam8 = lam08 / sy.sqrt(3)


#
def fsu3(i, j, k):
    """Calculate QCD structure constants
    by using permutations of know
    structure constants
    """
    sq3 = sy.sqrt(3) / 2 
    s12 = sy.Rational(1, 2)
    lf = [
        ((0, 1, 2), 1),
        ((1, 0, 2), -1),
        ((0, 3, 6), s12),
        ((3, 0, 6), -s12),
        ((1, 3, 5), s12),
        ((3, 1, 5), -s12),
        ((1, 4, 6), s12),
        ((4, 1, 6), -s12),
        ((2, 3, 4), s12),
        ((3, 2, 4), -s12),
        ((4, 0, 5), s12),
        ((0, 4, 5), -s12),
        ((5, 2, 6), s12),
        ((2, 5, 6), -s12),
        ((3, 4, 7), sq3),
        ((4, 3, 7), -sq3),
        ((5, 6, 7), sq3),
        ((6, 5, 7), -sq3),
    ]
    t = (i, j, k)
    #
    for tt in lf:
        if not set(t) == set(tt[0]):
            continue
        if Permutation([t]) == Permutation([tt[0]]):
            return tt[1]
        else:
            continue
    return 0


####
